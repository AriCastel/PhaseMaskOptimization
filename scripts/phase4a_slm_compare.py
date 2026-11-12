"""
Phase 4a: candidate masks on the real SLM, with a no-mask control.

SLM: 9.2 um pixels, mask radius 362 px (724 px across), 8-bit, phase range
[0, 2 pi) so the mask is Fresnel-cut (mod 2 pi). The back-focal-plane pupil is
3.30 mm in radius = 358.7 px, i.e. slightly smaller than the drawn mask.

Masks: no mask (control), best zoned_focus N = 3, best zoned_focus N = 4 + sweep,
best re-centred spherical (spherical + defocus), all from docs/phase2/results.json.

Device configurations (psfsim.slm_psf.DEVICES):
  ideal (continuous)                      reference = Phase 2 numbers
  SLM, pupil = mask (362 px)              pixels + 2 pi cut + 8 bit
  SLM, pupil 3.30 mm (358.7 px)           ... + real pupil size   <- "real SLM"
  SLM, pupil 3.30 mm + fringing 0.5 px    ... + pixel crosstalk (assumed)

Outputs in docs/phase4a/: results.json, key_depths.json, figures/*.png and the
8-bit SLM patterns (slm_patterns/*.png, 724 x 724, gray level = phase * 256 / 2pi).

Usage: python scripts/phase4a_slm_compare.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem, RadialPSF  # noqa: E402
from psfsim.crlb import Photometry  # noqa: E402
from psfsim.families import SphericalDefocus, ZonedFocus  # noqa: E402
from psfsim.slm_psf import DEVICES, REAL_PUPIL_PX, SLMPSF, range_summary_slm  # noqa: E402

OUT = ROOT / "docs" / "phase4a"
REAL = "SLM, pupil 3.30 mm (358.7 px)"
Z = np.arange(-4000.0, 4000.0 + 1, 50.0)
KEY_Z = [-4000.0, -3000.0, -2000.0, -1000.0, 0.0, 1000.0, 2000.0, 3000.0, 4000.0]


def candidate_masks():
    cfg = json.loads((ROOT / "docs" / "phase2" / "results.json").read_text())["configs"]
    return {
        "no mask (control)": (None, None),
        "zoned_focus N=3": ZonedFocus(3).mask(cfg["zoned_focus_3"]["x"]),
        "zoned_focus N=4 + sweep": ZonedFocus(4, sweep=True).mask(cfg["zoned_focus_4_sweep"]["x"]),
        "spherical + defocus (re-centred)": SphericalDefocus().mask(cfg["spherical+defocus"]["x"]),
    }


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "slm_patterns").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    system = OpticalSystem()
    masks = candidate_masks()
    budget = RadialPSF(system).budgets().power_ropt_nm
    print(f"SLM: mask radius 362 px, pupil {REAL_PUPIL_PX:.1f} px, 8 bit, 2pi cut; "
          f"camera {system.camera_pixel_object_nm:.1f} nm; z {Z[0]:.0f}..{Z[-1]:.0f} nm step 50", flush=True)

    sims = {name: SLMPSF(system, dev) for name, dev in DEVICES.items()}
    sweeps, summary = {}, {}
    for mname, (m, p) in masks.items():
        for dname, sim in sims.items():
            ph = sim.displayed_phase(m, p)
            sw = sim.sweep(Z, ph, keep_rows=(dname == REAL))
            sweeps[(mname, dname)] = sw
            s = range_summary_slm(sw, budget_nm=budget)
            for zz in (-4000.0, 4000.0):
                i = int(np.argmin(np.abs(Z - zz)))
                s[f"ee_at_{zz / 1000:+.0f}um"] = float(sw["encircled_energy"][i])
                s[f"strehl_at_{zz / 1000:+.0f}um"] = float(sw["strehl"][i])
                s[f"crlb_at_{zz / 1000:+.0f}um_nm"] = float(sw["crlb_nm"][i])
            summary.setdefault(mname, {})[dname] = s
            print(f"  {mname:34s} | {dname:38s} softmin {s['ee_softmin']:.4f} EEmin {s['ee_min']:.4f} "
                  f"Smin {s['strehl_min']:.4f} CRLBmax {s['crlb_max_nm']:.2f} nm ({time.time() - t0:.0f} s)",
                  flush=True)

    # ---- key depths on the real SLM -------------------------------------------
    sim = sims[REAL]
    rng = np.random.default_rng(2026)
    ph_nom = Photometry(camera_pixel_nm=system.camera_pixel_object_nm)
    key = {}
    for mname, (m, p) in masks.items():
        s = summary[mname][REAL]
        zs = KEY_Z + [s["z_crlb_worst_nm"], s["z_crlb_best_nm"]]
        labels = [f"{z / 1000:+.0f} µm" for z in KEY_Z] + [
            f"worst CRLB\n{s['z_crlb_worst_nm'] / 1000:+.2f} µm", f"best CRLB\n{s['z_crlb_best_nm'] / 1000:+.2f} µm"]
        phase = sim.displayed_phase(m, p)
        key[mname] = []
        for z, lab in zip(zs, labels):
            r = sim.plane(z, phase, keep_images=True)
            cam = ph_nom.n_photons * r["image_camera"] + ph_nom.background
            r["frame"] = rng.poisson(cam).astype(np.float32)
            r["expected"] = cam
            r["label"] = lab
            key[mname].append(r)

    # ---- tables -------------------------------------------------------------
    cols = [("ee_softmin", "soft-min EE"), ("ee_min", "EE min"), ("ee_mean", "EE mean"),
            ("strehl_min", "S min"), ("strehl_mean", "S mean"), ("crlb_max_nm", "worst CRLB"),
            ("crlb_mean_nm", "mean CRLB"), ("crlb_max_highbg_nm", "worst CRLB hi-bg"),
            ("ee50_max_core_nm", "EE50 max ±2"), ("feasible", "S≥0.03")]
    print("\n## Metrics over |z| <= 3 um, by SLM configuration\n")
    print("| mask | configuration | " + " | ".join(c[1] for c in cols) + " |")
    print("|---|---|" + "---|" * len(cols))
    for mname in masks:
        for dname in DEVICES:
            s = summary[mname][dname]
            v = []
            for k, _ in cols:
                x = s[k]
                v.append(("yes" if x else "no") if isinstance(x, bool) else
                         (f"{x:.0f}" if abs(x) >= 100 else f"{x:.3f}" if abs(x) < 10 else f"{x:.1f}"))
            print(f"| {mname} | {dname} | " + " | ".join(v) + " |")

    print("\n## Change from ideal to real SLM (pupil 3.30 mm), and with fringing\n")
    print("| mask | soft-min EE | EE min | S min | worst CRLB |")
    print("|---|---|---|---|---|")
    for mname in masks:
        a, b, c = (summary[mname][k] for k in ("ideal (continuous)", REAL, "SLM, pupil 3.30 mm + fringing 0.5 px"))
        f = lambda k: f"{100 * (b[k] / a[k] - 1):+.1f}% / {100 * (c[k] / a[k] - 1):+.1f}%"  # noqa: E731
        print(f"| {mname} | {f('ee_softmin')} | {f('ee_min')} | {f('strehl_min')} | {f('crlb_max_nm')} |")

    print("\n## Key depths on the real SLM (EE400 / peak Strehl / CRLB nm / EE50 nm)\n")
    hdr = [r["label"].replace("\n", " ") for r in key["no mask (control)"][:len(KEY_Z)]] + ["worst CRLB", "best CRLB"]
    print("| mask | " + " | ".join(hdr) + " |")
    print("|---|" + "---|" * len(hdr))
    for mname in masks:
        cells = []
        for r in key[mname]:
            extra = f" @{r['z_nm'] / 1000:+.2f}" if "CRLB" in r["label"] else ""
            cells.append(f"{r['encircled_energy']:.3f} / {r['strehl']:.3f} / {r['crlb_nm']:.1f} / "
                         f"{r['ee50_nm']:.0f}{extra}")
        print(f"| {mname} | " + " | ".join(cells) + " |")

    # ---- SLM patterns ---------------------------------------------------------
    from PIL import Image

    patterns = {}
    for mname, (m, p) in masks.items():
        lv = sim.slm_levels(m, p)
        patterns[mname] = lv
        fn = mname.split(" (")[0].replace(" ", "_").replace("=", "").replace("+", "plus")
        Image.fromarray(lv.astype(np.uint8), mode="L").save(OUT / "slm_patterns" / f"{fn}.png")
        np.save(OUT / "slm_patterns" / f"{fn}.npy", lv.astype(np.uint8))

    # ---- save -----------------------------------------------------------------
    js = {"slm": {"pitch_um": 9.2, "mask_radius_px": 362, "pupil_radius_px": REAL_PUPIL_PX, "bit_depth": 8,
                  "phase_range": "0..2pi (mod 2pi)", "subsample": sim.s},
          "photometry": ph_nom.__dict__, "masks": {k: repr(v) for k, v in masks.items()},
          "devices": {k: v.__dict__ for k, v in DEVICES.items()}, "summary": summary,
          "key_depths": {m: [{k: (float(v) if np.isscalar(v) and not isinstance(v, str) else v)
                              for k, v in r.items() if not isinstance(v, np.ndarray)} for r in rs]
                         for m, rs in key.items()},
          "runtime_s": time.time() - t0}
    (OUT / "results.json").write_text(json.dumps(js, indent=2, default=float))
    make_figures(masks, sweeps, summary, key, patterns, sim)
    print(f"\ndone in {time.time() - t0:.0f} s; wrote {OUT}")


def make_figures(masks, sweeps, summary, key, patterns, sim):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig_dir = OUT / "figures"
    names = list(masks)
    colors = {"no mask (control)": "#4a5568", "zoned_focus N=3": "#2f855a",
              "zoned_focus N=4 + sweep": "#2b6cb0", "spherical + defocus (re-centred)": "#c05621"}
    xf = sim.x_fine / 1000.0
    crop = np.abs(sim.x_fine) <= 2000
    ext = [xf[crop][0], xf[crop][-1], xf[crop][0], xf[crop][-1]]
    ncol = len(key[names[0]])

    # 1. fine PSF gallery (log scale, shared, relative to no-mask in-focus peak)
    fig, axes = plt.subplots(len(names), ncol, figsize=(2.0 * ncol, 2.25 * len(names) + 0.6))
    for i, mname in enumerate(names):
        for j, r in enumerate(key[mname]):
            ax = axes[i, j]
            img = r["image_fine"][np.ix_(crop, crop)] / sim.reference_peak
            im = ax.imshow(np.clip(img, 1e-4, None), extent=ext, origin="lower", cmap="magma",
                           norm=LogNorm(vmin=1e-3, vmax=1.0))
            circ = plt.Circle((0, 0), 0.4, fill=False, color="w", lw=0.6, ls="--")
            ax.add_patch(circ)
            ax.set_xticks([]), ax.set_yticks([])
            ax.set_title(r["label"], fontsize=7.5)
            ax.text(0.03, 0.03, f"EE {r['encircled_energy']:.2f}\nS {r['strehl']:.2f}\nσ {r['crlb_nm']:.1f} nm",
                    transform=ax.transAxes, fontsize=6, color="w", va="bottom")
            if j == 0:
                ax.set_ylabel(mname.replace(" (", "\n("), fontsize=8)
    fig.colorbar(im, ax=axes, shrink=0.6, label="intensity / no-mask in-focus peak")
    fig.suptitle("PSF on the real SLM (8 bit, 2π cut, 724 px, pupil 3.30 mm), fine grid ±2 µm, log scale; "
                 "dashed circle: r_opt = 400 nm; σ = lateral CRLB (3000 photons, 4.6 bg/px)", fontsize=10)
    fig.savefig(fig_dir / "psf_gallery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 2. camera frames (98.5 nm pixels, one Poisson realisation, per-row colour scale)
    fig, axes = plt.subplots(len(names), ncol, figsize=(2.0 * ncol, 2.25 * len(names) + 0.6))
    for i, mname in enumerate(names):
        vmax = max(r["frame"].max() for r in key[mname])
        for j, r in enumerate(key[mname]):
            ax = axes[i, j]
            ax.imshow(r["frame"], cmap="gray", vmin=0, vmax=vmax, origin="lower")
            ax.set_xticks([]), ax.set_yticks([])
            ax.set_title(r["label"], fontsize=7.5)
            ax.text(0.03, 0.03, f"peak {r['expected'].max():.0f} ph", transform=ax.transAxes, fontsize=6,
                    color="y", va="bottom")
            if j == 0:
                ax.set_ylabel(mname.replace(" (", "\n("), fontsize=8)
    fig.suptitle("Simulated camera frames, 21×21 px of 98.5 nm (2.07 µm), 3000 photons + 4.6 photons/px "
                 "background, one Poisson draw; colour scale shared within each row", fontsize=10)
    fig.savefig(fig_dir / "camera_gallery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 3. x-z sections on the real SLM
    fig, axes = plt.subplots(1, len(names), figsize=(4.0 * len(names), 4.6), sharey=True)
    for ax, mname in zip(axes, names):
        sw = sweeps[(mname, REAL)]
        img = sw["xz_fine"][:, crop] / sim.reference_peak
        im = ax.imshow(np.clip(img, 1e-4, None), extent=[xf[crop][0], xf[crop][-1], Z[0] / 1000, Z[-1] / 1000],
                       origin="lower", aspect="auto", cmap="magma", norm=LogNorm(vmin=1e-3, vmax=1.0))
        for y in (-3, 3):
            ax.axhline(y, color="w", lw=0.6, ls="--")
        for xx in (-0.4, 0.4):
            ax.axvline(xx, color="c", lw=0.5, ls=":")
        ax.set_title(mname, fontsize=9)
        ax.set_xlabel("x (µm)")
    axes[0].set_ylabel("z (µm)")
    fig.colorbar(im, ax=axes, shrink=0.8, label="intensity / no-mask in-focus peak")
    fig.suptitle("x–z sections on the real SLM (log scale); dashed: ±3 µm target range; dotted: r_opt", fontsize=10)
    fig.savefig(fig_dir / "xz_sections.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 4. metrics vs z, ideal vs SLM configurations
    styles = {"ideal (continuous)": (":", 1.2), "SLM, pupil = mask (362 px)": ("-.", 1.0),
              REAL: ("-", 1.8), "SLM, pupil 3.30 mm + fringing 0.5 px": ("--", 1.2)}
    fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)
    for mname in names:
        for dname, (ls, lw) in styles.items():
            sw = sweeps[(mname, dname)]
            kw = dict(color=colors[mname], ls=ls, lw=lw, label=f"{mname} — {dname}" if dname == REAL else None)
            axes[0].plot(Z / 1000, sw["encircled_energy"], **kw)
            axes[1].plot(Z / 1000, sw["strehl"], **kw)
            axes[2].plot(Z / 1000, sw["crlb_nm"], **kw)
    axes[0].axhspan(0.15, 0.19, color="#f6e05e", alpha=0.25)
    axes[1].axhline(0.03, color="k", lw=0.8, ls=":")
    axes[2].set_yscale("log")
    for ax, lab in zip(axes, ["EE400", "Strehl (peak)", "lateral CRLB (nm), 3000 ph, 4.6 bg/px"]):
        ax.set_ylabel(lab)
        ax.axvspan(-4, -3, color="0.94", zorder=0)
        ax.axvspan(3, 4, color="0.94", zorder=0)
        ax.grid(alpha=0.3)
    axes[2].set_xlabel("z (µm)")
    axes[0].legend(fontsize=7.5, loc="upper right")
    axes[0].set_title("Line style: solid = real SLM (pupil 3.30 mm); dotted = ideal continuous; dash-dot = SLM with "
                      "pupil = mask; dashed = real SLM + 0.5 px fringing", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig_dir / "metrics_vs_z.png", dpi=140)
    plt.close(fig)

    # 5. SLM patterns
    fig, axes = plt.subplots(1, len(names) - 1, figsize=(4.2 * (len(names) - 1), 4.4))
    for ax, mname in zip(axes, names[1:]):
        im = ax.imshow(patterns[mname], cmap="twilight", vmin=0, vmax=255, origin="lower")
        ax.add_patch(plt.Circle((361.5, 361.5), REAL_PUPIL_PX, fill=False, color="w", lw=0.8, ls="--"))
        ax.set_title(mname, fontsize=9)
        ax.set_xticks([0, 362, 723]), ax.set_yticks([0, 362, 723])
    fig.colorbar(im, ax=axes, shrink=0.8, label="8-bit level (256 = 2π)")
    fig.suptitle("8-bit SLM patterns (724 × 724 px, 9.2 µm); dashed: 3.30 mm pupil", fontsize=10)
    fig.savefig(fig_dir / "slm_patterns.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
