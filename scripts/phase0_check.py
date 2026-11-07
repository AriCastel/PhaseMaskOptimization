"""
Phase 0: fast evaluator and budget check.

1. Cross-check the radial (1D Debye) model against the 2D PSFSimulator for
   no mask, the 2 rad RMS spherical mask and the zoned +-1.5 pi mask over
   z = -8..+8 um (2D at 100 nm steps; radial at 50 nm).
2. Check numerically that Integral S_axial dz and Integral P(r_opt) dz are the
   same for every radially symmetric mask, and equal the analytic budgets.
3. Report, for each mask: in-range fraction eta over |z| <= 3 um, EE floor and
   flatness, mean Strehl, and the budget-limited flat-EE ceiling.

Usage:
    python scripts/phase0_check.py               # full run (~2-3 min, 2D dominates)
    python scripts/phase0_check.py --quick       # 2D at 21 planes only
    python scripts/phase0_check.py --apodization uniform

Writes docs/phase0/results.json and docs/phase0/figures/*.png.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem, PSFSimulator, RadialPSF  # noqa: E402
from psfsim.metrics import plane_metrics, range_summary  # noqa: E402

R_OPT = 400.0
TOTAL_R = 10_000.0
Z_RANGE = (-3000.0, 3000.0)
SPH_2RAD = 2.0 / (2.0 * np.pi)  # 2 rad RMS in waves (Noll-normalised Z11)

# name -> (mask, params)
CROSSCHECK_MASKS = {
    "none": (None, None),
    "spherical 2 rad RMS": ("spherical", {"strength_waves": SPH_2RAD}),
    "zoned +-1.5pi": ("zoned_zernike", {"zone_coefficients": ({4: 0.75}, {4: -0.75})}),
}
BUDGET_MASKS = {
    **CROSSCHECK_MASKS,
    "zoned +-0.5pi": ("zoned_zernike", {"zone_coefficients": ({4: 0.25}, {4: -0.25})}),
    "zoned +-pi": ("zoned_zernike", {"zone_coefficients": ({4: 0.5}, {4: -0.5})}),
    "annular 0/pi/0": ("annular", None),
    "defocus_sweep 2 um": ("defocus_sweep", {"depth_range_nm": 2000.0}),
    "random m=0 Zernikes": ("zernike", {"coefficients": {11: 0.4, 22: -0.3, 37: 0.2}}),
}
METRICS = ("strehl", "strehl_axial", "encircled_energy", "ee50_nm", "power_weighted_m2", "power_ropt")


def run_2d(sim: PSFSimulator, z: np.ndarray, mask, params, ref_power: float) -> dict:
    out = {k: [] for k in METRICS}
    px = sim.system.sim_pixel_nm
    c = sim.system.grid_size // 2
    for r in sim.iter_stack(z, mask=mask, mask_params=params, slm_pixelation=False):
        m = plane_metrics(r.psf, px, R_OPT, TOTAL_R, ref_power)
        out["strehl"].append(r.strehl)
        out["strehl_axial"].append(r.psf[c, c] / r._reference_peak)
        out["encircled_energy"].append(m.encircled_energy)
        out["ee50_nm"].append(r.encircled_energy_radius_nm)
        out["power_weighted_m2"].append(m.power_weighted_m2)
        out["power_ropt"].append(m.power_ropt)
    return {k: np.asarray(v) for k, v in out.items()}


def agreement(rad_z, rad_v, z2, v2, z_lim=None):
    """Max |2D - radial| relative to the radial curve's maximum (on the 2D planes)."""
    rv = np.interp(z2, rad_z, rad_v)
    sel = np.ones_like(z2, bool) if z_lim is None else (np.abs(z2) <= z_lim)
    scale = np.nanmax(np.abs(rv[sel]))
    return float(np.nanmax(np.abs(v2[sel] - rv[sel])) / scale)


def pointwise(rad_z, rad_v, z2, v2, z_lim=None):
    """Max pointwise relative error |2D/radial - 1| on the 2D planes."""
    rv = np.interp(z2, rad_z, rad_v)
    sel = np.ones_like(z2, bool) if z_lim is None else (np.abs(z2) <= z_lim)
    return float(np.nanmax(np.abs(v2[sel] / rv[sel] - 1.0)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--apodization", default="uniform")
    ap.add_argument("--out", default=str(ROOT / "docs" / "phase0"))
    args = ap.parse_args()
    out_dir = Path(args.out)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    # Continuous phase (no 8-bit quantisation, no pixelation) so the two
    # models see exactly the same pupil; the SLM is Phase 4.
    system = OpticalSystem(apodization=args.apodization, slm_bit_depth=None)
    rad = RadialPSF(system, r_opt_nm=R_OPT, total_radius_nm=TOTAL_R)
    budgets = rad.budgets()
    z_rad = np.arange(-8000.0, 8000.0 + 1, 50.0)
    z_wide = np.arange(-60000.0, 60000.0 + 1, 25.0)
    z_2d = np.linspace(-5000, 5000, 21) if args.quick else np.arange(-8000.0, 8000.0 + 1, 100.0)

    print(f"system: NA {system.na}, {system.wavelength_nm} nm, n {system.n_medium}, "
          f"apodization {system.apodization}; 2D grid {system.grid_size}, "
          f"{system.sim_pixel_nm:.1f} nm/px")
    print(f"analytic budgets: Int S_axial dz = {budgets.strehl_axial_nm:.1f} nm "
          f"(uniform-g formula {budgets.extra['uniform_g_strehl_nm']:.1f}, paraxial "
          f"{budgets.extra['paraxial_strehl_nm']:.1f}); Int P(r_opt) dz = {budgets.power_ropt_nm:.1f} nm")

    # ---------------- radial sweeps + budgets for every mask --------------
    results = {"system": {"na": system.na, "wavelength_nm": system.wavelength_nm,
                          "n_medium": system.n_medium, "apodization": system.apodization,
                          "r_opt_nm": R_OPT, "total_radius_nm": TOTAL_R},
               "budgets": {"strehl_axial_nm": budgets.strehl_axial_nm,
                           "power_ropt_nm": budgets.power_ropt_nm, **budgets.extra},
               "masks": {}}
    stacks = {}
    t0 = time.time()
    for name, (mask, params) in BUDGET_MASKS.items():
        st = rad.simulate(z_rad, mask, params, label=name)
        stacks[name] = st
        phase = rad.total_phase(mask, params)
        s_wide = rad.axial_strehl(z_wide, phase)
        p_wide = rad.intensity(z_wide, phase)[:, : rad._i_opt + 1]
        rr = rad.r_nm[: rad._i_opt + 1]
        p_wide = simpson(p_wide * 2 * np.pi * rr, x=rr, axis=1)
        ee = range_summary(z_rad, st["encircled_energy"], Z_RANGE)
        pr = range_summary(z_rad, st["power_ropt"], Z_RANGE, total=budgets.power_ropt_nm)
        sa = range_summary(z_rad, st["strehl_axial"], Z_RANGE, total=budgets.strehl_axial_nm)
        sp = range_summary(z_rad, st["strehl"], Z_RANGE)
        i0 = int(np.argmin(np.abs(z_rad)))
        results["masks"][name] = {
            "mask": mask, "params": repr(params),
            "int_strehl_axial_8um": float(np.trapezoid(st["strehl_axial"], z_rad)),
            "int_strehl_axial_60um": float(np.trapezoid(s_wide, z_wide)),
            "int_strehl_peak_8um": float(np.trapezoid(st["strehl"], z_rad)),
            "int_power_ropt_8um": float(np.trapezoid(st["power_ropt"], z_rad)),
            "int_power_ropt_60um": float(np.trapezoid(p_wide, z_wide)),
            "int_ee_8um": float(np.trapezoid(st["encircled_energy"], z_rad)),
            "strehl_z0": float(st["strehl"][i0]), "ee_z0": float(st["encircled_energy"][i0]),
            "ee_range": ee, "power_ropt_range": pr, "strehl_axial_range": sa,
            "strehl_peak_range": sp,
        }
    print(f"radial sweeps: {time.time() - t0:.1f} s for {len(BUDGET_MASKS)} masks "
          f"x {z_rad.size} planes (+ wide-range budget checks)")

    print("\n## Budget check (nm)\n")
    print("| mask | Int S_ax ±8 | Int S_ax ±60 | Int P(r_opt) ±8 | Int P(r_opt) ±60 | Int S_peak ±8 | Int EE ±8 |")
    print("|---|---|---|---|---|---|---|")
    for name, m in results["masks"].items():
        print(f"| {name} | {m['int_strehl_axial_8um']:.1f} | {m['int_strehl_axial_60um']:.1f} | "
              f"{m['int_power_ropt_8um']:.1f} | {m['int_power_ropt_60um']:.1f} | "
              f"{m['int_strehl_peak_8um']:.1f} | {m['int_ee_8um']:.1f} |")
    print(f"| analytic (infinite z) | {budgets.strehl_axial_nm:.1f} | {budgets.strehl_axial_nm:.1f} | "
          f"{budgets.power_ropt_nm:.1f} | {budgets.power_ropt_nm:.1f} | not invariant | not invariant |")

    print("\n## In-range summary, |z| <= 3 um (radial model)\n")
    print("| mask | S(0) | EE(0) | eta_P (vs budget) | eta_EE (vs ±8) | EE min | EE mean | flatness | mean S_peak |")
    print("|---|---|---|---|---|---|---|---|---|")
    for name, m in results["masks"].items():
        e, p, s = m["ee_range"], m["power_ropt_range"], m["strehl_peak_range"]
        print(f"| {name} | {m['strehl_z0']:.3f} | {m['ee_z0']:.3f} | {p['eta']:.3f} | {e['eta']:.3f} | "
              f"{e['min']:.3f} | {e['mean']:.3f} | {e['flatness']:.2f} | {s['mean']:.3f} |")
    ceil_p = budgets.power_ropt_nm / (Z_RANGE[1] - Z_RANGE[0])
    print(f"\nflat ceiling over ±3 um: P(r_opt) <= {ceil_p:.3f}; "
          f"over ±4 um: {budgets.power_ropt_nm / 8000:.3f} (fractions of total light)")
    results["ceiling_power_ropt_3um"] = ceil_p
    results["ceiling_power_ropt_4um"] = budgets.power_ropt_nm / 8000

    # ---------------- 2D cross-check ---------------------------------------
    sim = PSFSimulator(system)
    ref2d = plane_metrics(sim.simulate(z_nm=0.0).psf, system.sim_pixel_nm, R_OPT, TOTAL_R).power_ropt
    cross = {}
    t0 = time.time()
    for name, (mask, params) in CROSSCHECK_MASKS.items():
        cross[name] = run_2d(sim, z_2d, mask, params, ref2d)
    print(f"\n2D sweeps: {time.time() - t0:.1f} s for {len(CROSSCHECK_MASKS)} masks x {z_2d.size} planes")

    print("\n## Radial vs 2D: max |diff| / max(curve)  [all z | |z|<=3 um]\n")
    print("| mask | " + " | ".join(METRICS) + " |")
    print("|---|" + "---|" * len(METRICS))
    results["crosscheck"] = {"z_2d_nm": z_2d.tolist(), "masks": {}}
    for name in CROSSCHECK_MASKS:
        st, c2 = stacks[name], cross[name]
        row, entry = [], {}
        for k in METRICS:
            a_all = agreement(z_rad, st[k], z_2d, c2[k])
            a_in = agreement(z_rad, st[k], z_2d, c2[k], 3000)
            pw_in = pointwise(z_rad, st[k], z_2d, c2[k], 3000)
            entry[k] = {"maxdiff_rel_all": a_all, "maxdiff_rel_3um": a_in,
                        "pointwise_rel_3um": pw_in, "values_2d": c2[k].tolist()}
            row.append(f"{100*a_all:.2f}% / {100*a_in:.2f}%")
        results["crosscheck"]["masks"][name] = entry
        print(f"| {name} | " + " | ".join(row) + " |")
        b2s = float(np.trapezoid(c2["strehl_axial"], z_2d))
        b2p = float(np.trapezoid(c2["power_ropt"], z_2d))
        entry["int_strehl_axial_2d"] = b2s
        entry["int_power_ropt_2d"] = b2p
    print("\n2D budget integrals over the 2D z grid (trapezoid):")
    for name in CROSSCHECK_MASKS:
        e = results["crosscheck"]["masks"][name]
        print(f"  {name:24s} Int S_axial = {e['int_strehl_axial_2d']:.1f} nm, "
              f"Int P(r_opt) = {e['int_power_ropt_2d']:.1f} nm")

    # ---------------- apodization sensitivity ------------------------------
    results["apodization"] = {}
    for apo in ("emission", "uniform", "illumination"):
        r_ = RadialPSF(system.with_(apodization=apo), r_opt_nm=R_OPT, total_radius_nm=TOTAL_R)
        b = r_.budgets()
        st0 = r_.simulate([0.0])
        results["apodization"][apo] = {"strehl_axial_nm": b.strehl_axial_nm,
                                       "power_ropt_nm": b.power_ropt_nm,
                                       "ee_z0_nomask": float(st0["encircled_energy"][0]),
                                       "ee50_z0_nomask": float(st0["ee50_nm"][0])}
    print("\n## Apodization sensitivity (no mask)\n")
    print("| apodization | Int S_ax dz (nm) | Int P(r_opt) dz (nm) | EE400(0) | EE50(0) nm |")
    print("|---|---|---|---|---|")
    for apo, v in results["apodization"].items():
        print(f"| {apo} | {v['strehl_axial_nm']:.0f} | {v['power_ropt_nm']:.0f} | "
              f"{v['ee_z0_nomask']:.3f} | {v['ee50_z0_nomask']:.0f} |")

    (out_dir / "results.json").write_text(json.dumps(results, indent=2))
    make_figures(out_dir / "figures", z_rad, stacks, z_2d, cross, budgets, ceil_p, results["masks"])
    print(f"\nwrote {out_dir / 'results.json'} and figures/")


def make_figures(fig_dir, z_rad, stacks, z_2d, cross, budgets, ceil_p, results_masks):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = ["#2b6cb0", "#c05621", "#2f855a", "#805ad5", "#b83280", "#718096", "#d69e2e", "#319795"]
    zr = z_rad / 1000

    # 1. cross-check
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    panels = [("strehl", "Strehl (peak)"), ("encircled_energy", "EE(400 nm)"),
              ("ee50_nm", "EE50 radius (nm)"), ("power_weighted_m2", "PW-M2 (nm²)")]
    for ax, (k, lab) in zip(axes.ravel(), panels):
        for col, name in zip(colors, cross):
            ax.plot(zr, stacks[name][k], color=col, lw=1.6, label=f"{name} (radial)")
            ax.plot(z_2d / 1000, cross[name][k], "o", ms=3.2, mfc="none", color=col,
                    label=f"{name} (2D)")
        ax.set_ylabel(lab)
        ax.axvspan(-3, 3, color="0.92", zorder=0)
        ax.grid(alpha=0.3)
        if k in ("ee50_nm", "power_weighted_m2"):
            ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("z (µm)")
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("Phase 0 cross-check: radial model (lines) vs 2D PSFSimulator (circles)")
    fig.tight_layout()
    fig.savefig(fig_dir / "crosscheck.png", dpi=150)
    plt.close(fig)

    # 2. budgets: where each mask puts the (shared) P(r_opt) budget
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    show = ["none", "spherical 2 rad RMS", "zoned +-1.5pi", "zoned +-pi"]
    for col, name in zip(colors, show):
        ax.plot(zr, stacks[name]["power_ropt"], color=col, lw=1.6, label=name)
    ax.axhline(ceil_p, color="k", ls="--", lw=1)
    ax.text(-7.8, ceil_p + 0.01, f"flat ceiling over ±3 µm = {ceil_p:.3f}", fontsize=8)
    ax.axvspan(-3, 3, color="0.92", zorder=0)
    ax.set_ylim(0, 0.5)
    ax.set_xlabel("z (µm)")
    ax.set_ylabel("P(r_opt = 400 nm), fraction of light")
    ax.set_title("P(r_opt) vs depth (no-mask peak 0.88 clipped)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    names = list(results_masks)
    inside = [results_masks[n]["power_ropt_range"]["integral"] for n in names]
    total = [results_masks[n]["int_power_ropt_60um"] for n in names]
    y = np.arange(len(names))
    ax.barh(y, inside, color="#2b6cb0", label="inside |z| ≤ 3 µm")
    ax.barh(y, np.subtract(total, inside), left=inside, color="#a0aec0",
            label="outside (to ±60 µm)")
    ax.axvline(budgets.power_ropt_nm, color="k", ls="--", lw=1)
    ax.text(budgets.power_ropt_nm, -0.75, f"analytic budget {budgets.power_ropt_nm:.0f} nm ",
            fontsize=8, va="bottom", ha="right")
    for yi, (a, t) in enumerate(zip(inside, total)):
        ax.text(10, yi, f"η = {a / budgets.power_ropt_nm:.2f}", va="center", fontsize=8, color="white")
    ax.set_yticks(y, names, fontsize=8)
    ax.set_xlabel("∫ P(r_opt) dz (nm)")
    ax.set_title("Same total for every mask; masks only move it in z", fontsize=10)
    ax.set_ylim(len(names) - 0.5, -1.0)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    fig.tight_layout()
    fig.savefig(fig_dir / "budgets.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
