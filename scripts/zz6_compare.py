"""
PSF comparison: the best zoned radial-Zernike mask (free Z(n, 0), n <= 6 per zone)
against the current zoned_focus N = 4 masks, the optimised spherical mask and no mask.

Masks
  no mask (control)
  spherical + defocus (re-centred)       optimised spherical aberration, Phase 2
  zoned_focus N=4 (defocus + piston)     Phase 2
  zoned_focus N=4 + sweep                Phase 2 (bench candidate)
  zoned_zernike6 N=<3|4>                 best of docs/zoned_zernike6/emitted/results.json, i.e. optimised
                                         with the 400 nm energy taken as a fraction of all emitted light
The mask optimised with the agreed EE400 (fraction of the light inside 10 um) is shown
only in figures/loophole.png: it wins that score by throwing light out of the window.

SLM configurations (psfsim.slm_psf, pixel-level model of Phase 4a)
  ideal (continuous)                          = the radial model of the search
  SLM, mask on the pupil (358.7 px)           pixels + 2 pi cut + 8 bit, mask drawn at the
                                              measured pupil radius (Phase 4a recommendation)  <- headline
  SLM, mask 362 px on a 358.7 px pupil        the 0.9 % size mismatch of Phase 4a
  SLM, mask on the pupil + fringing 0.5 px    pixel crosstalk (assumed width)

Writes docs/zoned_zernike6/: compare.json, figures/*.png, slm_patterns/.

Usage: python scripts/zz6_compare.py [--workers 4]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem, RadialPSF  # noqa: E402
from psfsim.crlb import Photometry  # noqa: E402
from psfsim.families import SphericalDefocus, ZonedFocus, ZonedZernikeRadial  # noqa: E402
from psfsim.objective import Evaluator, Targets, soft_min, trapezoid_weights  # noqa: E402
from psfsim.slm_psf import REAL_PUPIL_PX, SLMDevice, SLMPSF, range_summary_slm  # noqa: E402

OUT = ROOT / "docs" / "zoned_zernike6"
Z = np.arange(-4000.0, 4000.0 + 1, 50.0)
KEY_Z = [-3000.0, -2000.0, -1000.0, 0.0, 1000.0, 2000.0, 3000.0]

HEAD = "SLM, mask on the pupil (358.7 px)"
DEVICES = {
    "ideal (continuous)": SLMDevice(ideal=True),
    HEAD: SLMDevice(mask_radius_px=REAL_PUPIL_PX, pupil_radius_px=REAL_PUPIL_PX),
    "SLM, mask 362 px on a 358.7 px pupil": SLMDevice(pupil_radius_px=REAL_PUPIL_PX),
    "SLM, mask on the pupil + fringing 0.5 px": SLMDevice(mask_radius_px=REAL_PUPIL_PX,
                                                         pupil_radius_px=REAL_PUPIL_PX, fringe_sigma_px=0.5),
}

# One colour per mask, fixed (validated categorical order: blue, orange, aqua, yellow; control in gray).
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
NEW, ZF4, ZF4S, SPH, NOMASK = "new", "zf4", "zf4s", "sph", "none"
COLOR = {NEW: "#2a78d6", ZF4: "#eb6834", ZF4S: "#1baf7a", SPH: "#eda100", NOMASK: "#8a8984"}


def best_new(denominator="emitted"):
    cfg = json.loads((OUT / denominator / "results.json").read_text())["configs"]
    cands = [k for k in ("zoned_zernike6_3", "zoned_zernike6_4") if cfg[k]["feasible"]]
    k = max(cands, key=lambda k: cfg[k]["score"])
    return k, cfg[k]


def candidate_masks():
    p2 = json.loads((ROOT / "docs" / "phase2" / "results.json").read_text())["configs"]
    k, e = best_new()
    n = int(k.rsplit("_", 1)[1])
    return {
        NOMASK: ("no mask (control)", (None, None)),
        SPH: ("spherical + defocus (optimised)", SphericalDefocus().mask(p2["spherical+defocus"]["x"])),
        ZF4: ("zoned defocus N=4 (defocus + piston)", ZonedFocus(4).mask(p2["zoned_focus_4"]["x"])),
        ZF4S: ("zoned defocus N=4 + sweep", ZonedFocus(4, sweep=True).mask(p2["zoned_focus_4_sweep"]["x"])),
        NEW: (f"zoned Zernike ≤ 6, N={n} (new)", ZonedZernikeRadial(n, 6).mask(e["x"])),
    }


def _sweep(args):
    key, dname = args
    _, (m, p) = candidate_masks()[key]
    sim = SLMPSF(OpticalSystem(), DEVICES[dname])
    return key, dname, sim.sweep(Z, sim.displayed_phase(m, p), keep_rows=(dname == HEAD))


def summarise(sw, budget):
    """range_summary_slm plus the emitted-light normalisation of the 400 nm energy."""
    s = range_summary_slm(sw, budget_nm=budget)
    sel = np.abs(sw["z_nm"]) <= 3000 + 1e-9
    p = sw["power_ropt"][sel]
    s.update({"p_softmin": soft_min(p, -10.0, trapezoid_weights(sw["z_nm"][sel])), "p_min": float(p.min()),
              "window_fraction_min": float(sw["total_power"][sel].min())})
    return s


def slm_indicators(m, p, radius_px=REAL_PUPIL_PX):
    """Phase jumps at the zone edges (size after the 2 pi cut, waves), the steepest slope
    inside a zone (waves per SLM pixel) and the narrowest zone (SLM pixels)."""
    rad = RadialPSF(n_pupil=1 << 16)
    order = np.argsort(rad.rho)
    rho, ph = rad.rho[order], rad.mask_phase(m, p)[order]
    edges = np.array(p.get("radii", ())) if p else np.array([])
    d = np.diff(ph)
    at_edge = np.zeros(d.size, bool)
    for e in edges:
        at_edge[np.searchsorted(rho, e) - 1] = True
    jumps = [float(abs((v + 0.5) % 1.0 - 0.5)) for v in d[at_edge]]
    slope = np.abs(d[~at_edge] / np.diff(rho)[~at_edge]).max() / radius_px if d.size else 0.0
    rings = np.diff(np.concatenate([[0.0], edges, [1.0]])) * radius_px
    return {"jumps_waves": jumps, "max_slope_waves_per_px": float(slope),
            "narrowest_zone_px": float(rings.min()), "peak_to_valley_waves": float(np.ptp(ph))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "slm_patterns").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    system = OpticalSystem()
    masks = candidate_masks()
    budget = RadialPSF(system).budgets().power_ropt_nm

    sweeps, summary = {}, {}
    with Pool(args.workers) as pool:
        for key, dname, sw in pool.imap_unordered(_sweep, [(k, d) for k in masks for d in DEVICES]):
            sweeps[(key, dname)] = sw
            summary.setdefault(key, {})[dname] = summarise(sw, budget)
            print(f"  {masks[key][0]:38s} | {dname:42s} softmin P {summary[key][dname]['p_softmin']:.4f} "
                  f"({time.time() - t0:.0f} s)", flush=True)

    # ---- key depths, headline configuration ------------------------------------
    sim = SLMPSF(system, DEVICES[HEAD])
    rng = np.random.default_rng(2026)
    ph_nom = Photometry(camera_pixel_nm=system.camera_pixel_object_nm)
    key = {}
    for k, (label, (m, p)) in masks.items():
        s = summary[k][HEAD]
        phase = sim.displayed_phase(m, p)
        key[k] = []
        for z, lab in zip(KEY_Z + [s["z_crlb_worst_nm"]],
                          [f"{z / 1000:+.0f} µm" for z in KEY_Z] + [f"worst CRLB: {s['z_crlb_worst_nm'] / 1000:+.2f} µm"]):
            r = sim.plane(z, phase, keep_images=True)
            cam = ph_nom.n_photons * r["image_camera"] + ph_nom.background
            r.update(frame=rng.poisson(cam).astype(np.float32), expected=cam, label=lab)
            key[k].append(r)

    # ---- tables ---------------------------------------------------------------
    cols = [("p_softmin", "soft-min P400 (emitted)"), ("p_min", "P400 min"), ("ee_softmin", "soft-min EE400 (window)"),
            ("ee_min", "EE400 min"), ("window_fraction_min", "light inside 10 µm, min"), ("strehl_min", "Strehl min"), ("crlb_max_nm", "worst CRLB"), ("crlb_mean_nm", "mean CRLB"),
            ("crlb_max_highbg_nm", "worst CRLB, high bg"), ("ee50_max_core_nm", "EE50 max ±2 µm")]

    def fmt(x):
        return f"{x:.0f}" if abs(x) >= 100 else f"{x:.3f}" if abs(x) < 1 else f"{x:.1f}"

    print("\n## Metrics over |z| <= 3 um, by SLM configuration\n")
    print("| mask | configuration | " + " | ".join(c[1] for c in cols) + " |")
    print("|---|---|" + "---|" * len(cols))
    for k, (label, _) in masks.items():
        for dname in DEVICES:
            print(f"| {label} | {dname} | " + " | ".join(fmt(summary[k][dname][c]) for c, _ in cols) + " |")

    print("\n## Key depths, " + HEAD + " (P400 / peak Strehl / CRLB nm / EE50 nm)\n")
    print("| mask | " + " | ".join(f"{z / 1000:+.0f} µm" for z in KEY_Z) + " | worst CRLB |")
    print("|---|" + "---|" * (len(KEY_Z) + 1))
    for k, (label, _) in masks.items():
        cells = [f"{r['power_ropt']:.3f} / {r['strehl']:.3f} / {r['crlb_nm']:.1f} / {r['ee50_nm']:.0f}"
                 + (f" @ {r['z_nm'] / 1000:+.2f}" if "worst" in r["label"] else "") for r in key[k]]
        print(f"| {label} | " + " | ".join(cells) + " |")

    ind = {k: slm_indicators(m, p) for k, (_, (m, p)) in masks.items() if m}
    print("\n## SLM indicators (mask radius 358.7 px)\n")
    print("| mask | phase jumps after the 2π cut (waves) | steepest slope (waves / px) | narrowest zone (px) |")
    print("|---|---|---|---|")
    for k, v in ind.items():
        print(f"| {masks[k][0]} | {', '.join(f'{j:.2f}' for j in v['jumps_waves']) or '–'} | "
              f"{v['max_slope_waves_per_px']:.3f} | {v['narrowest_zone_px']:.0f} |")

    # ---- SLM pattern of the new mask (drawn at the pupil radius, 724 x 724 frame) --------
    from PIL import Image

    patterns = {}
    for k, (label, (m, p)) in masks.items():
        lv = sim.slm_levels(m, p).astype(np.uint8)
        patterns[k] = np.pad(lv, (724 - lv.shape[0]) // 2)
    name = best_new()[0]
    Image.fromarray(patterns[NEW], mode="L").save(OUT / "slm_patterns" / f"{name}_r358.7.png")
    np.save(OUT / "slm_patterns" / f"{name}_r358.7.npy", patterns[NEW])

    js = {"masks": {k: {"label": lab, "mask": repr(mp)} for k, (lab, mp) in masks.items()},
          "headline": HEAD, "devices": {k: v.__dict__ for k, v in DEVICES.items()}, "summary": summary,
          "slm_indicators": ind, "photometry": ph_nom.__dict__,
          "key_depths": {k: [{a: (float(b) if np.isscalar(b) and not isinstance(b, str) else b)
                              for a, b in r.items() if not isinstance(b, np.ndarray)} for r in rs]
                         for k, rs in key.items()},
          "runtime_s": time.time() - t0}
    (OUT / "compare.json").write_text(json.dumps(js, indent=2, default=float))
    np.savez_compressed(OUT / "compare_curves.npz", z_nm=Z, **{
        f"{k}__{d}__{m}": sweeps[(k, d)][m] for k in masks for d in DEVICES
        for m in ("encircled_energy", "power_ropt", "total_power", "strehl", "crlb_nm", "crlb_highbg_nm", "ee50_nm")})
    make_figures(masks, sweeps, summary, key, patterns, sim)
    print(f"\ndone in {time.time() - t0:.0f} s; wrote {OUT}")


def _style():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK, "xtick.color": INK2,
        "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.titlelocation": "left", "lines.linewidth": 2.0,
        "lines.solid_capstyle": "round", "legend.frameon": False, "xtick.major.size": 0, "ytick.major.size": 0,
    })
    return plt


def make_figures(masks, sweeps, summary, key, patterns, sim):
    plt = _style()
    from matplotlib.colors import LogNorm

    fig_dir = OUT / "figures"
    names = list(masks)
    lab = {k: v[0] for k, v in masks.items()}
    zum = Z / 1000.0

    def line_kw(k):
        return dict(color=COLOR[k], lw=2.6 if k == NEW else 1.8, ls="--" if k == NOMASK else "-",
                    zorder=5 if k == NEW else 3, label=lab[k])

    # 1. metrics vs z (headline configuration): three small multiples on a shared z axis
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.5), sharex=True)
    panels = [("power_ropt", "P400: fraction of all emitted light within r = 400 nm"),
              ("strehl", "Peak Strehl (log scale)"),
              ("crlb_nm", "Lateral CRLB, nm (log scale; 3000 photons, 4.6 photons/px background)")]
    for ax, (m, title) in zip(axes, panels):
        for k in names:
            ax.plot(zum, sweeps[(k, HEAD)][m], **line_kw(k))
        ax.set_title(title)
        ax.axvspan(-4, -3, color="#f0efec", zorder=0, lw=0)
        ax.axvspan(3, 4, color="#f0efec", zorder=0, lw=0)
    axes[0].set_ylim(0, 0.45)
    axes[0].annotate("no mask peaks at 0.88", xy=(0.12, 0.445), xytext=(0.5, 0.40), fontsize=8, color=INK2,
                     arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    axes[1].set_yscale("log")
    axes[1].set_ylim(3e-3, 1.2)
    axes[1].axhline(0.03, color=INK2, lw=0.8)
    axes[1].text(-3.95, 0.033, "constraint: Strehl ≥ 0.03", fontsize=8, color=INK2, va="bottom")
    axes[2].set_yscale("log")
    axes[2].set_ylim(1.5, 150)
    axes[2].set_yticks([2, 5, 10, 20, 50, 100], labels=["2", "5", "10", "20", "50", "100"])
    axes[2].set_xlabel("z (µm); shaded: outside the ±3 µm design range")
    axes[0].legend(ncol=2, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    fig.suptitle("Depth response of the five PSFs (" + HEAD + ")", x=0.012, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(fig_dir / "metrics_vs_z.png", dpi=150)
    plt.close(fig)

    # 2. worst case within the range: one horizontal bar chart per metric
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.3))
    for ax, (m, title, f) in zip(axes, [("p_softmin", "Soft-min P400 over ±3 µm (higher is better)", "{:.3f}"),
                                        ("p_min", "Lowest P400 over ±3 µm (higher is better)", "{:.3f}"),
                                        ("crlb_max_nm", "Worst lateral CRLB over ±3 µm, nm (lower is better)", "{:.1f}")]):
        order = names[::-1]
        vals = [summary[k][HEAD][m] for k in order]
        ax.barh(range(len(order)), vals, height=0.55, color=[COLOR[k] for k in order])
        for i, v in enumerate(vals):
            ax.text(v, i, " " + f.format(v), va="center", fontsize=8.5, color=INK)
        ax.set_yticks(range(len(order)), [lab[k] for k in order] if ax is axes[0] else [""] * len(order))
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", visible=False)
        ax.set_xlim(0, max(vals) * 1.18)
    fig.tight_layout()
    fig.savefig(fig_dir / "worst_case.png", dpi=150)
    plt.close(fig)

    # 3. PSF gallery (fine grid, log scale, relative to the no-mask in-focus peak)
    xf = sim.x_fine / 1000.0
    crop = np.abs(sim.x_fine) <= 2000
    ext = [xf[crop][0], xf[crop][-1], xf[crop][0], xf[crop][-1]]
    ncol = len(key[names[0]])
    fig, axes = plt.subplots(len(names), ncol, figsize=(1.95 * ncol + 1.2, 2.1 * len(names) + 0.7))
    for i, k in enumerate(names):
        for j, r in enumerate(key[k]):
            ax = axes[i, j]
            img = r["image_fine"][np.ix_(crop, crop)] / sim.reference_peak
            im = ax.imshow(np.clip(img, 1e-4, None), extent=ext, origin="lower", cmap="magma",
                           norm=LogNorm(vmin=1e-3, vmax=1.0))
            ax.add_patch(plt.Circle((0, 0), 0.4, fill=False, color="w", lw=0.6, ls="--"))
            ax.set_xticks([]), ax.set_yticks([])
            ax.grid(False)
            ax.set_title(r["label"], fontsize=7.5, fontweight="normal", loc="center")
            ax.text(0.03, 0.03, f"P {r['power_ropt']:.2f}\nS {r['strehl']:.2f}\nσ {r['crlb_nm']:.1f} nm",
                    transform=ax.transAxes, fontsize=6.5, color="w", va="bottom")
            if j == 0:
                ax.set_ylabel(lab[k].replace(" (", "\n(").replace(", N", ",\nN"), fontsize=8, color=INK)
    fig.colorbar(im, ax=axes, shrink=0.5, pad=0.015, label="intensity / no-mask in-focus peak")
    fig.suptitle("PSF at each depth (" + HEAD + "; 4 × 4 µm, log scale). Dashed circle: r = 400 nm. "
                 "P = P400, S = peak Strehl, σ = lateral CRLB", fontsize=10, x=0.01, ha="left")
    fig.savefig(fig_dir / "psf_gallery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 4. simulated camera frames (one Poisson draw; one grey scale for all masked rows)
    vmax = max(r["frame"].max() for k in names if k != NOMASK for r in key[k])
    fig, axes = plt.subplots(len(names), ncol, figsize=(1.95 * ncol + 1.2, 2.1 * len(names) + 0.7))
    for i, k in enumerate(names):
        for j, r in enumerate(key[k]):
            ax = axes[i, j]
            im = ax.imshow(r["frame"], cmap="gray", vmin=0, vmax=vmax, origin="lower")
            ax.set_xticks([]), ax.set_yticks([])
            ax.grid(False)
            ax.set_title(r["label"], fontsize=7.5, fontweight="normal", loc="center")
            ax.text(0.03, 0.03, f"peak {r['expected'].max():.0f} ph", transform=ax.transAxes, fontsize=6.5,
                    color="#f5d76e", va="bottom")
            if j == 0:
                ax.set_ylabel(lab[k].replace(" (", "\n(").replace(", N", ",\nN"), fontsize=8, color=INK)
    fig.colorbar(im, ax=axes, shrink=0.5, pad=0.015, label="photons per pixel (clipped at the top of the scale)")
    fig.suptitle("Simulated camera frames: 21 × 21 px of 98.5 nm, 3000 photons + 4.6 photons/px background, "
                 "one Poisson draw, one grey scale for all panels", fontsize=10, x=0.01, ha="left")
    fig.savefig(fig_dir / "camera_gallery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 5. x-z sections
    fig, axes = plt.subplots(1, len(names), figsize=(3.1 * len(names) + 1, 4.8), sharey=True)
    for ax, k in zip(axes, names):
        img = sweeps[(k, HEAD)]["xz_fine"][:, crop] / sim.reference_peak
        im = ax.imshow(np.clip(img, 1e-4, None), extent=[xf[crop][0], xf[crop][-1], zum[0], zum[-1]],
                       origin="lower", aspect="auto", cmap="magma", norm=LogNorm(vmin=1e-3, vmax=1.0))
        for y in (-3, 3):
            ax.axhline(y, color="w", lw=0.6, ls="--")
        ax.grid(False)
        ax.set_title(lab[k].replace(" (", "\n("), fontsize=8.5)
        ax.set_xlabel("x (µm)")
    axes[0].set_ylabel("z (µm)")
    fig.colorbar(im, ax=axes, shrink=0.8, pad=0.015, label="intensity / no-mask in-focus peak")
    fig.suptitle("x–z sections through the axis (log scale); dashed: the ±3 µm design range", fontsize=10,
                 x=0.01, ha="left")
    fig.savefig(fig_dir / "xz_sections.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 6. phase profiles and SLM patterns of the masks
    rad = RadialPSF(n_pupil=4096)
    shown = [k for k in names if k != NOMASK]
    fig = plt.figure(figsize=(12.5, 6.6))
    gs = fig.add_gridspec(2, len(shown), height_ratios=[1.0, 1.25])
    ax = fig.add_subplot(gs[0, :])
    order = np.argsort(rad.rho)
    for k in shown:
        ph = rad.mask_phase(*masks[k][1])[order]
        ax.plot(rad.rho[order], ph - ph[0], **line_kw(k))
    ax.set_xlabel("normalised pupil radius ρ")
    ax.set_title("Mask phase before the 2π cut (waves; centre = 0)")
    ax.legend(ncol=4, fontsize=8.5, loc="upper left")
    for j, k in enumerate(shown):
        a = fig.add_subplot(gs[1, j])
        im = a.imshow(patterns[k], cmap="twilight", vmin=0, vmax=255, origin="lower")
        a.set_title(lab[k], fontsize=8.5)
        a.set_xticks([0, 362, 723]), a.set_yticks([0, 362, 723])
        a.grid(False)
    fig.colorbar(im, ax=fig.axes[1:], shrink=0.75, pad=0.015, label="8-bit level (256 = 2π)")
    fig.suptitle("Phase profiles and 8-bit SLM patterns (mask radius 358.7 px, 724 × 724 px frame)", fontsize=11,
                 x=0.01, ha="left", fontweight="bold")
    fig.savefig(fig_dir / "phase_and_patterns.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 7. search: best score of every optimiser run, per family
    js = json.loads((OUT / "emitted" / "results.json").read_text())["configs"]
    fams = [("zoned_focus_3", "zoned defocus\nN=3", ZF4), ("zoned_focus_3_sweep", "zoned defocus\nN=3 + sweep", ZF4S),
            ("zoned_zernike6_3", "zoned Zernike ≤ 6\nN=3", NEW), ("zoned_focus_4", "zoned defocus\nN=4", ZF4),
            ("zoned_focus_4_sweep", "zoned defocus\nN=4 + sweep", ZF4S), ("zoned_zernike6_4", "zoned Zernike ≤ 6\nN=4", NEW)]
    mk = {"cold": "o", "cold-long": "s", "warm": "^", "polish": "D", "phase2": "*"}
    lo_y = max(js[f]["search_score"] for f, _, _ in fams) - 0.035
    fig, ax = plt.subplots(figsize=(10, 4.8))
    rng = np.random.default_rng(0)
    for i, (f, _, ck) in enumerate(fams):
        for r in js[f]["runs"]:
            if r["score"] < lo_y:
                continue
            src = r["source"].split(":")[0]
            ax.scatter(i + rng.uniform(-0.22, 0.22), r["score"], s=26, marker=mk[src], color=COLOR[ck], alpha=0.75,
                       edgecolor=SURFACE, linewidth=0.6, zorder=3)
        ax.hlines(js[f]["search_score"], i - 0.32, i + 0.32, color=INK, lw=1.2, zorder=4)
        ax.text(i, js[f]["search_score"] + 0.0012, f"{js[f]['search_score']:.4f}", ha="center", fontsize=8.5)
        if "phase2_score" in js[f]:
            ax.scatter(i - 0.38, js[f]["phase2_score"], s=70, marker="*", color=INK2, zorder=5)
    for src, label in [("cold", "random start, Phase 2 budget"), ("cold-long", "random start, 3× population, 12 000 evaluations"),
                       ("warm", "started from a zoned-defocus optimum or the window-objective optimum"), ("polish", "polish of the best runs"),
                       ("phase2", "Phase 2 optimum of that family, same score")]:
        ax.scatter([], [], s=70 if src == "phase2" else 26, marker=mk[src], color=INK2, label=label)
    ax.axvline(2.5, color=GRID, lw=1)
    ax.set_xticks(range(len(fams)), [f[1] for f in fams])
    top = max(js[f]["search_score"] for f, _, _ in fams)
    ax.set_ylim(top - 0.035, top + 0.006)
    ax.set_ylabel("soft-min P400 over ±3 µm (Strehl ≥ 0.03)")
    ax.grid(axis="x", visible=False)
    ax.set_title("Best score of every CMA-ES run (search evaluator); bar and number: best of the family.\n"
                 "Runs below the axis range are not shown.")
    ax.legend(fontsize=7.5, loc="lower left", ncol=2)
    fig.tight_layout()
    fig.savefig(fig_dir / "search_runs.png", dpi=150)
    plt.close(fig)

    # 8. the loophole: EE400 (window) against P400 (emitted) for the window-optimised mask
    ev = Evaluator(RadialPSF(n_pupil=2048, dr_nm=10.0), Targets(dz_nm=25.0, z_range_nm=(-4000.0, 4000.0)))
    kw, ew = best_new("window")
    nw = int(kw.rsplit("_", 1)[1])
    trio = [(ZF4S, lab[ZF4S], masks[ZF4S][1], "-"),
            (NEW, lab[NEW].replace(" (new)", ", optimised on P400 (emitted)"), masks[NEW][1], "-"),
            ("loop", f"zoned Zernike ≤ 6, N={nw}, optimised on EE400 (window)", ZonedZernikeRadial(nw, 6).mask(ew["x"]), "-")]
    col = {**COLOR, "loop": "#e34948"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), sharex=True)
    for k, label, (m, p), ls in trio:
        r = ev.evaluate_mask(m, p, keep_curves=True)
        zz = r["curve_z_nm"] / 1000
        for ax, c in zip(axes, ("curve_encircled_energy", "curve_power_ropt", "curve_total_power")):
            ax.plot(zz, r[c], color=col[k], lw=2.0, ls=ls, label=label)
    for ax, t in zip(axes, ["EE400: energy within 400 nm / light inside 10 µm\n(the agreed score)",
                            "P400: energy within 400 nm / all emitted light",
                            "Fraction of the emitted light inside 10 µm\n(the EE400 denominator)"]):
        ax.set_title(t, fontsize=9)
        ax.set_xlabel("z (µm)")
        ax.axvspan(-4, -3, color="#f0efec", zorder=0, lw=0)
        ax.axvspan(3, 4, color="#f0efec", zorder=0, lw=0)
    axes[0].set_ylim(0, 0.34), axes[1].set_ylim(0, 0.34), axes[2].set_ylim(0, 1.05)
    axes[0].legend(fontsize=7.5, loc="lower center")
    fig.tight_layout()
    fig.savefig(fig_dir / "loophole.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
