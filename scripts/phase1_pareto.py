"""
Phase 1: objective, targets and a small Pareto study.

Two mask families (both radially symmetric, both contain the baselines):
  A  "zernike_m0": Noll-normalised Z4 (defocus), Z11, Z22, Z37 (RMS waves)
  B  "zoned2":     two zones, edge rho1, plain defocus c_in / c_out (peak waves,
                   the zoned_zernike convention) and a piston on the outer zone

Study 1 (trade-off, as in the plan): NSGA-II on
    max EE_min(|z|<=3),  max mean Strehl(|z|<=3),  max eta
Study 2 (thresholds): NSGA-II on
    max softmin EE,  max Strehl_min(|z|<=3),  min EE50_max(|z|<=2)
  Its front answers "best attainable EE floor for given (S_min, R_max)" exactly,
  which is what the advisor needs to fix the thresholds.

Every evaluation (params + all metrics) is logged to docs/phase1/evals.csv.

Usage:  python scripts/phase1_pareto.py [--pop 100 --gen 80]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pymoo.algorithms.moo.nsga2 import NSGA2  # noqa: E402
from pymoo.core.problem import Problem  # noqa: E402
from pymoo.optimize import minimize  # noqa: E402
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting  # noqa: E402

from psfsim import OpticalSystem, RadialPSF  # noqa: E402
from psfsim.objective import Evaluator, Targets  # noqa: E402

OUT = ROOT / "docs" / "phase1"
SPH_2RAD = 1.0 / np.pi

FAMILIES = {
    "zernike_m0": {
        "names": ["z4", "z11", "z22", "z37"],
        "lo": [-1.2, -1.0, -0.6, -0.4],
        "hi": [1.2, 1.0, 0.6, 0.4],
    },
    "zoned2": {
        "names": ["rho1", "c_in", "c_out", "piston"],
        "lo": [0.3, -1.5, -1.5, 0.0],
        "hi": [0.95, 1.5, 1.5, 1.0],
    },
}

STUDY_OBJECTIVES = {
    "tradeoff": [("ee_min", -1), ("strehl_mean", -1), ("eta", -1)],
    "thresholds": [("ee_softmin", -1), ("strehl_min", -1), ("ee50_max_core_nm", +1)],
}

SUMMARY_KEYS = ["ee_softmin", "ee_min", "ee_mean", "ee_flatness", "ee_cv", "eta",
                "strehl_min", "strehl_mean", "ee50_max_core_nm", "pwm2_mean_nm2",
                "violation", "feasible", "score"]


def family_mask(family: str, x) -> tuple:
    if family == "zernike_m0":
        z4, z11, z22, z37 = x
        return "zernike", {"coefficients": {4: z4, 11: z11, 22: z22, 37: z37}}
    rho1, c_in, c_out, pist = x
    return "zoned_zernike", {"radii": (float(rho1),),
                             "zone_coefficients": ({4: c_in}, {4: c_out, 1: pist})}


class Log:
    def __init__(self):
        self.rows = []

    def add(self, source, family, x, m):
        row = {"source": source, "family": family}
        row.update({f"x{i}": float(v) for i, v in enumerate(x)})
        row.update({k: m[k] for k in SUMMARY_KEYS})
        if "crlb_max_nm" in m:
            row["crlb_max_nm"] = m["crlb_max_nm"]
        self.rows.append(row)

    def write(self, path):
        keys = sorted({k for r in self.rows for k in r}, key=lambda k: (k not in ("source", "family"), k))
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(self.rows)


class MaskProblem(Problem):
    def __init__(self, ev, family, objectives, log, tag):
        spec = FAMILIES[family]
        super().__init__(n_var=len(spec["lo"]), n_obj=len(objectives),
                         xl=np.array(spec["lo"]), xu=np.array(spec["hi"]))
        self.ev, self.family, self.objectives, self.log, self.tag = ev, family, objectives, log, tag

    def _evaluate(self, X, out, *args, **kwargs):
        F = np.empty((X.shape[0], self.n_obj))
        for i, x in enumerate(X):
            m = self.ev.evaluate_mask(*family_mask(self.family, x))
            self.log.add(self.tag, self.family, x, m)
            F[i] = [sgn * min(m[k], 1e6) for k, sgn in self.objectives]
        out["F"] = F


def run_nsga(ev, family, study, log, pop, gen, seed):
    prob = MaskProblem(ev, family, STUDY_OBJECTIVES[study], log, f"nsga_{study}")
    res = minimize(prob, NSGA2(pop_size=pop), ("n_gen", gen), seed=seed, verbose=False)
    return res


def baselines(ev):
    """Named reference masks, plus spherical re-centred / re-tuned with defocus."""
    base = {
        "none": (None, None),
        "spherical 2 rad RMS": ("spherical", {"strength_waves": SPH_2RAD}),
        "zoned ±0.5π": ("zoned_zernike", {"zone_coefficients": ({4: 0.25}, {4: -0.25})}),
        "zoned ±π": ("zoned_zernike", {"zone_coefficients": ({4: 0.5}, {4: -0.5})}),
        "zoned ±1.5π": ("zoned_zernike", {"zone_coefficients": ({4: 0.75}, {4: -0.75})}),
        "annular 0/π/0": ("annular", None),
        "defocus_sweep 2 um": ("defocus_sweep", {"depth_range_nm": 2000.0}),
    }
    # spherical 2 rad RMS + best defocus offset (Z4), by grid
    z4s = np.linspace(-1.0, 1.0, 201)
    sc = [ev.evaluate_mask("zernike", {"coefficients": {4: a, 11: SPH_2RAD}})["ee_softmin"] for a in z4s]
    a_best = float(z4s[int(np.argmax(sc))])
    base["spherical 2 rad + defocus (re-centred)"] = ("zernike", {"coefficients": {4: a_best, 11: SPH_2RAD}})
    # spherical strength and defocus both tuned, by grid
    best = (-1, None)
    for s11 in np.linspace(0.05, 0.8, 31):
        for a in np.linspace(-1.0, 1.0, 81):
            v = ev.evaluate_mask("zernike", {"coefficients": {4: a, 11: s11}})["ee_softmin"]
            if v > best[0]:
                best = (v, (a, s11))
    a, s11 = best[1]
    base["spherical + defocus (both tuned)"] = ("zernike", {"coefficients": {4: a, 11: s11}})
    return base


def front_indices(F):
    return NonDominatedSorting().do(F, only_non_dominated_front=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", type=int, default=100)
    ap.add_argument("--gen", type=int, default=80)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    (OUT / "figures").mkdir(parents=True, exist_ok=True)

    system = OpticalSystem()
    ev = Evaluator(RadialPSF(system), Targets())
    log = Log()
    t0 = time.time()
    info = ev.describe()
    print(f"targets: {info['targets']}")
    print(f"P(r_opt) budget {info['budget_power_ropt_nm']:.0f} nm -> flat ceiling {info['flat_ceiling']:.3f}")

    # ---- baselines --------------------------------------------------------
    base = baselines(ev)
    base_res = {}
    for name, (mask, params) in base.items():
        m = ev.evaluate_mask(mask, params, with_crlb=True, keep_curves=True)
        base_res[name] = m
        log.add("baseline:" + name, "baseline", [], m)
    print(f"baselines done ({time.time() - t0:.0f} s)")

    # ---- NSGA-II studies ---------------------------------------------------
    fronts = {}
    for study in STUDY_OBJECTIVES:
        for fam in FAMILIES:
            res = run_nsga(ev, fam, study, log, args.pop, args.gen, args.seed)
            fronts[(study, fam)] = res
            print(f"NSGA-II {study:10s} {fam:10s}: {len(res.F)} front points "
                  f"({time.time() - t0:.0f} s, {ev.n_evals} evals)")

    # ---- re-evaluate the trade-off fronts with curves + CRLB ---------------
    picks = {}
    tradeoff_pts = []
    for fam in FAMILIES:
        res = fronts[("tradeoff", fam)]
        for x in res.X:
            m = ev.evaluate_mask(*family_mask(fam, x), with_crlb=True)
            log.add("front_tradeoff", fam, x, m)
            tradeoff_pts.append((fam, x, m))
    # union front
    F = np.array([[-m["ee_min"], -m["strehl_mean"], -m["eta"]] for _, _, m in tradeoff_pts])
    union = front_indices(F)
    tradeoff_pts = [tradeoff_pts[i] for i in union]

    # ---- thresholds: attainable floor vs (S_min, R_max) -------------------
    thr_pts = []
    for fam in FAMILIES:
        res = fronts[("thresholds", fam)]
        for x in res.X:
            m = ev.evaluate_mask(*family_mask(fam, x))
            thr_pts.append((fam, x, m))
    # all logged NSGA points also count (they are real evaluations)
    pool = [(r["family"], [r[f"x{i}"] for i in range(4)], r) for r in log.rows
            if r["source"].startswith("nsga")] + thr_pts
    s_grid = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    r_grid = [1000, 1500, 2000, 2500, 3000, 4000, np.inf]
    table = np.full((len(s_grid), len(r_grid)), np.nan)
    table_min = np.full_like(table, np.nan)
    best_at = {}
    for i, smin in enumerate(s_grid):
        for j, rmax in enumerate(r_grid):
            ok = [p for p in pool if p[2]["strehl_min"] >= smin and p[2]["ee50_max_core_nm"] <= rmax]
            if ok:
                b = max(ok, key=lambda p: p[2]["ee_softmin"])
                table[i, j] = b[2]["ee_softmin"]
                table_min[i, j] = b[2]["ee_min"]
                best_at[(smin, rmax)] = b
    print("\n## Best attainable softmin EE400 over |z|<=3 um (rows S_min, cols R_max nm)\n")
    print("| S_min \\ R_max | " + " | ".join("∞" if np.isinf(r) else f"{r:.0f}" for r in r_grid) + " |")
    print("|---|" + "---|" * len(r_grid))
    for i, smin in enumerate(s_grid):
        print(f"| {smin:.2f} | " + " | ".join("–" if np.isnan(v) else f"{v:.3f}" for v in table[i]) + " |")

    # ---- picks --------------------------------------------------------------
    def pick(name, fam, x):
        m = ev.evaluate_mask(*family_mask(fam, x), with_crlb=True, keep_curves=True)
        picks[name] = {"family": fam, "x": [float(v) for v in x], "metrics": m}

    fam, x, _ = max(tradeoff_pts, key=lambda p: p[2]["ee_min"])
    pick("max EE_min (trade-off front)", fam, x)
    fam, x, _ = max(thr_pts, key=lambda p: p[2]["ee_softmin"])
    pick("max softmin EE (unconstrained)", fam, x)
    tg = ev.targets
    for fam_name in FAMILIES:
        cands = [p for p in pool if p[0] == fam_name and p[2]["strehl_min"] >= tg.s_min
                 and p[2]["ee50_max_core_nm"] <= tg.r_max_nm]
        if cands:
            b = max(cands, key=lambda p: p[2]["ee_softmin"])
            pick(f"best feasible ({tg.s_min}, {tg.r_max_nm:.0f} nm), {fam_name}", fam_name, b[1])

    # ---- tables -------------------------------------------------------------
    cols = ["ee_softmin", "ee_min", "ee_mean", "eta", "ee_flatness", "strehl_min", "strehl_mean",
            "ee50_max_core_nm", "crlb_max_nm", "feasible"]
    print("\n## Baselines and picks (|z|<=3 um; feasibility at provisional targets)\n")
    print("| mask | " + " | ".join(cols) + " |")
    print("|---|" + "---|" * len(cols))

    def fmt(v):
        if isinstance(v, (bool, np.bool_)):
            return "yes" if v else "no"
        return f"{v:.0f}" if abs(v) >= 100 else f"{v:.3f}"

    for name, m in list(base_res.items()) + [(k, v["metrics"]) for k, v in picks.items()]:
        print(f"| {name} | " + " | ".join(fmt(m[c]) for c in cols) + " |")
    for name, v in picks.items():
        print(f"  {name}: {v['family']} x = {np.round(v['x'], 3).tolist()}")
    for name in ("spherical 2 rad + defocus (re-centred)", "spherical + defocus (both tuned)"):
        print(f"  {name}: {base[name][1]}")

    # ---- does EE ranking agree with the CRLB objective? --------------------
    ee_v = np.array([m["ee_softmin"] for _, _, m in tradeoff_pts])
    ee_m = np.array([m["ee_min"] for _, _, m in tradeoff_pts])
    cr_v = np.array([m["crlb_max_nm"] for _, _, m in tradeoff_pts])
    from scipy.stats import spearmanr
    rho_s = spearmanr(ee_v, -cr_v).statistic
    rho_m = spearmanr(ee_m, -cr_v).statistic
    print(f"\nSpearman(softmin EE, -worst CRLB) over {len(cr_v)} trade-off front points: {rho_s:.2f}; "
          f"(EE_min, -worst CRLB): {rho_m:.2f}")

    # ---- save -----------------------------------------------------------------
    log.write(OUT / "evals.csv")
    js = {"describe": info, "s_grid": s_grid, "r_grid": [None if np.isinf(r) else r for r in r_grid],
          "best_softmin_table": table.tolist(), "best_eemin_table": table_min.tolist(),
          "spearman_softmin_vs_crlb": rho_s, "spearman_eemin_vs_crlb": rho_m,
          "baselines": {k: {kk: (float(vv) if np.isscalar(vv) else None) for kk, vv in m.items()
                            if not kk.startswith("curve")} for k, m in base_res.items()},
          "baseline_params": {k: repr(v) for k, v in base.items()},
          "picks": {k: {"family": v["family"], "x": v["x"],
                        "metrics": {kk: (float(vv) if np.isscalar(vv) else None)
                                    for kk, vv in v["metrics"].items() if not kk.startswith("curve")}}
                    for k, v in picks.items()},
          "n_evals": ev.n_evals}
    (OUT / "results.json").write_text(json.dumps(js, indent=2, default=float))
    make_figures(ev, base_res, picks, tradeoff_pts, thr_pts, table, s_grid, r_grid, info)
    print(f"\n{ev.n_evals} evaluations in {time.time() - t0:.0f} s; wrote {OUT}")


def make_figures(ev, base_res, picks, tradeoff_pts, thr_pts, table, s_grid, r_grid, info):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = OUT / "figures"
    ceil = info["flat_ceiling"]
    marks = {"zernike_m0": "o", "zoned2": "^"}

    # 1. trade-off front
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (xk, xl) in zip(axes, [("strehl_mean", "mean Strehl over |z| ≤ 3 µm"),
                                   ("eta", "η (in-range share of P(r_opt) budget)")]):
        for fam, mk in marks.items():
            pts = [m for f, _, m in tradeoff_pts if f == fam]
            sc = ax.scatter([m[xk] for m in pts], [m["ee_min"] for m in pts], c=[m["eta" if xk != "eta" else "strehl_mean"] for m in pts],
                            marker=mk, s=22, cmap="viridis", label=f"front: {fam}")
        for name, m in base_res.items():
            ax.scatter(m[xk], m["ee_min"], marker="x", color="#c05621", s=40)
            ax.annotate(name, (m[xk], m["ee_min"]), fontsize=6.5, xytext=(3, 3),
                        textcoords="offset points", color="#c05621")
        ax.axhline(ceil, color="k", ls="--", lw=1)
        ax.text(ax.get_xlim()[0], ceil * 1.02, f" flat ceiling {ceil:.3f}", fontsize=8)
        ax.set_xlabel(xl)
        ax.set_ylabel("EE400 min over |z| ≤ 3 µm")
        ax.grid(alpha=0.3)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("η" if xk != "eta" else "mean Strehl")
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("Phase 1 trade-off: NSGA-II fronts (two mask families) and baselines (×)")
    fig.tight_layout()
    fig.savefig(fig_dir / "pareto_tradeoff.png", dpi=150)
    plt.close(fig)

    # 2. threshold map
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    im = ax.imshow(table, origin="lower", cmap="viridis", aspect="auto", vmin=0)
    ax.set_xticks(range(len(r_grid)), ["∞" if np.isinf(r) else f"{r:.0f}" for r in r_grid])
    ax.set_yticks(range(len(s_grid)), [f"{s:.2f}" for s in s_grid])
    for i in range(len(s_grid)):
        for j in range(len(r_grid)):
            v = table[i, j]
            ax.text(j, i, "–" if np.isnan(v) else f"{v:.3f}", ha="center", va="center",
                    fontsize=8, color="w" if (np.isnan(v) or v < np.nanmax(table) * 0.6) else "k")
    ax.set_xlabel("R_max: EE50 cap within |z| ≤ 2 µm (nm)")
    ax.set_ylabel("S_min: Strehl floor within |z| ≤ 3 µm")
    ax.set_title(f"Best attainable soft-min EE400 over ±3 µm (ceiling {ceil:.3f})", fontsize=10)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(fig_dir / "thresholds.png", dpi=150)
    plt.close(fig)

    # 3. curves
    show = {k: base_res[k] for k in ("spherical 2 rad RMS", "spherical 2 rad + defocus (re-centred)",
                                     "spherical + defocus (both tuned)", "zoned ±π")}
    show.update({k: v["metrics"] for k, v in picks.items()})
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for col, (name, m) in zip(colors, show.items()):
        z = m["curve_z_nm"] / 1000
        ls = "--" if name in base_res else "-"
        axes[0, 0].plot(z, m["curve_encircled_energy"], ls, color=col, label=name)
        axes[0, 1].plot(z, m["curve_strehl"], ls, color=col)
        axes[1, 0].plot(z, m["curve_ee50_nm"], ls, color=col)
        axes[1, 1].plot(z, m["curve_crlb_nm"], ls, color=col)
    axes[0, 0].axhline(ceil, color="k", lw=0.8, ls=":")
    axes[0, 1].axhline(ev.targets.s_min, color="k", lw=0.8, ls=":")
    axes[1, 0].axhline(ev.targets.r_max_nm, color="k", lw=0.8, ls=":")
    axes[1, 0].axvspan(-3, -2, color="0.93", zorder=0)
    axes[1, 0].axvspan(2, 3, color="0.93", zorder=0)
    labels = ["EE400", "Strehl (peak)", "EE50 (nm)", "lateral CRLB (nm, placeholder photometry)"]
    for ax, lab in zip(axes.ravel(), labels):
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[1, 1].set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("z (µm)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Baselines (dashed) and Phase 1 picks (solid); dotted: ceiling / provisional targets")
    fig.tight_layout()
    fig.savefig(fig_dir / "curves.png", dpi=150)
    plt.close(fig)

    # 4. EE vs CRLB
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for fam, mk in marks.items():
        pts = [m for f, _, m in tradeoff_pts if f == fam]
        ax.scatter([m["ee_softmin"] for m in pts], [m["crlb_max_nm"] for m in pts], marker=mk, s=18,
                   label=fam)
    for name, m in base_res.items():
        ax.scatter(m["ee_softmin"], m["crlb_max_nm"], marker="x", color="#c05621")
    ax.set_yscale("log")
    ax.set_xlabel("soft-min EE400 over ±3 µm")
    ax.set_ylabel("worst lateral CRLB over ±3 µm (nm)")
    ax.set_title("Primary vs alternative objective (trade-off front points)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "ee_vs_crlb.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
