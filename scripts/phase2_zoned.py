"""
Phase 2: zoned family, re-parameterised (`zoned_focus`), against re-optimised baselines.

Objective (agreed after Phase 1): maximise the soft-min of EE400 over |z| <= 3 um
subject to peak Strehl >= 0.03 over the same range; EE50 is only reported.

Searches (all on the fast radial evaluator: 1024 pupil rings, dr 20 nm, dz 50 nm):
  * dense grid for zoned_focus N = 2 (edge, two depths, piston), then local CMA-ES
    polishing of the best grid points, with and without the sweep term;
  * CMA-ES, 20 random restarts, for zoned_focus N = 2..5, with and without sweeps;
  * the same CMA-ES protocol (10-20 restarts) for every baseline family:
    spherical + defocus (re-centred spherical), Z4/Z11/Z22/Z37 (Phase 1 best family),
    the current two-zone zoned_zernike, binary 0/pi annular masks with 1-3 edges,
    and the coherent layer cake (equal dc zones, pistons only) for N = 2..5;
  * Optuna TPE (Bayesian optimisation) as an independent cross-check on
    zoned_focus N = 2, 3 and 3 + sweep.
Every evaluation is logged to docs/phase2/evals.csv.gz. The best mask of every
configuration is then re-scored on the reference evaluator (2048 rings, dr 10 nm,
dz 12.5 nm) with the lateral CRLB in the nominal and high-background regimes.

Usage:  python scripts/phase2_zoned.py [--restarts 20] [--workers 4] [--quick]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import gzip  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem, RadialPSF  # noqa: E402
from psfsim.crlb import HIGH_BACKGROUND, LateralCRLB, Photometry  # noqa: E402
from psfsim.families import (AnnularBinary, LayerCake, SphericalDefocus, ZernikeM0,  # noqa: E402
                             ZonedFocus, ZonedZernike2)
from psfsim.objective import Evaluator, Targets  # noqa: E402

OUT = ROOT / "docs" / "phase2"
LOG_KEYS = ["score", "ee_softmin", "ee_min", "ee_mean", "eta", "ee_flatness", "strehl_min",
            "strehl_mean", "ee50_max_core_nm", "pwm2_mean_nm2", "feasible"]

_EV = None


def _init_worker():
    global _EV
    _EV = Evaluator(RadialPSF(n_pupil=1024, dr_nm=20.0), Targets())


def build_family(spec):
    kind, *args = spec
    return {"sph": SphericalDefocus, "zm0": ZernikeM0, "zz2": ZonedZernike2, "ann": AnnularBinary,
            "cake": LayerCake, "zf": ZonedFocus}[kind](*args)


def _row(source, fam, run, x, m):
    r = {"source": source, "family": fam.name, "run": run,
         "x": json.dumps([round(float(v), 6) for v in x])}
    r.update({k: (bool(m[k]) if k == "feasible" else round(float(m[k]), 6)) for k in LOG_KEYS})
    return r


# --------------------------------------------------------------------------- #
# Tasks (run in workers)
# --------------------------------------------------------------------------- #
def task_cma(args):
    import cma

    spec, seed, budget, x0_u, sigma0, source = args
    fam = build_family(spec)
    rng = np.random.default_rng(seed)
    u0 = rng.uniform(0.05, 0.95, fam.dim) if x0_u is None else np.clip(np.asarray(x0_u), 0, 1)
    span = fam.hi - fam.lo
    if fam.dim == 1:  # CMA-ES needs >= 2 dimensions: dense 1-D scan instead
        rows, best = [], (-np.inf, None)
        for u in np.linspace(0, 1, budget) if seed == 0 else []:
            x = fam.lo + u * span
            m = _EV.evaluate_mask(*fam.mask(x))
            rows.append(_row("scan", fam, seed, x, m))
            if m["score"] > best[0]:
                best = (m["score"], x)
        return {"spec": spec, "source": "cma", "seed": seed, "best_score": best[0],
                "best_x": best[1] if best[1] is not None else fam.lo.copy(), "rows": rows}
    es = cma.CMAEvolutionStrategy(u0, sigma0, {"bounds": [0, 1], "seed": int(seed) + 1, "maxfevals": budget,
                                               "verbose": -9, "tolfun": 1e-6, "tolx": 1e-6})
    rows, best = [], (-np.inf, None)
    while not es.stop():
        U = es.ask()
        F = []
        for u in U:
            x = fam.lo + np.clip(u, 0, 1) * span
            m = _EV.evaluate_mask(*fam.mask(x))
            rows.append(_row(source, fam, seed, x, m))
            F.append(-m["score"])
            if m["score"] > best[0]:
                best = (m["score"], x)
        es.tell(U, F)
    return {"spec": spec, "source": source, "seed": seed, "best_score": best[0], "best_x": best[1], "rows": rows}


def task_grid(args):
    spec, xs = args
    fam = build_family(spec)
    rows = []
    for x in xs:
        m = _EV.evaluate_mask(*fam.mask(x))
        rows.append(_row("grid", fam, 0, x, m))
    return {"spec": spec, "source": "grid", "rows": rows}


def task_optuna(args):
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    spec, n_trials, seed = args
    fam = build_family(spec)
    rows = []

    def objective(trial):
        x = np.array([trial.suggest_float(n, lo, hi) for n, lo, hi in zip(fam.names, fam.lo, fam.hi)])
        m = _EV.evaluate_mask(*fam.mask(x))
        rows.append(_row("optuna_tpe", fam, seed, x, m))
        return m["score"]

    study = optuna.create_study(direction="maximize",
                                sampler=optuna.samplers.TPESampler(seed=seed, multivariate=True))
    study.optimize(objective, n_trials=n_trials)
    bx = np.array([study.best_params[n] for n in fam.names])
    return {"spec": spec, "source": "optuna_tpe", "seed": seed, "best_score": study.best_value,
            "best_x": bx, "rows": rows}


# --------------------------------------------------------------------------- #
def budget_for(dim):
    return int(np.clip(250 * dim, 1000, 4000))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--quick", action="store_true", help="tiny run for testing the pipeline")
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from saved results")
    args = ap.parse_args()
    if args.figures_only:
        js = json.loads((OUT / "results.json").read_text())
        z = np.load(OUT / "best_curves.npz")
        curves = {}
        for key in z.files:
            k, kk = key.split("__")
            curves.setdefault(k, {})[kk] = z[key]
        make_figures(js["configs"], curves, js["order"], Evaluator(RadialPSF(n_pupil=256), Targets()))
        return
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    R = 2 if args.quick else args.restarts
    bscale = 0.1 if args.quick else 1.0
    t0 = time.time()

    baseline_specs = [("sph",), ("zm0",), ("zz2",), ("ann", 1), ("ann", 2), ("ann", 3),
                      ("cake", 2), ("cake", 3), ("cake", 4), ("cake", 5)]
    zoned_specs = [("zf", n, sw) for n in (2, 3, 4, 5) for sw in (False, True)]

    tasks = []
    for spec in baseline_specs:
        fam = build_family(spec)
        nr = R if spec[0] == "zm0" else max(2, R // 2)
        tasks += [(task_cma, (spec, s, int(budget_for(fam.dim) * bscale), None, 0.3, "cma")) for s in range(nr)]
    for spec in zoned_specs:
        fam = build_family(spec)
        tasks += [(task_cma, (spec, s, int(budget_for(fam.dim) * bscale), None, 0.3, "cma")) for s in range(R)]
    # Optuna cross-check
    n_trials = 150 if args.quick else 2000
    tasks += [(task_optuna, (spec, n_trials, 7)) for spec in [("zf", 2, False), ("zf", 3, False), ("zf", 3, True)]]
    # dense grid for N = 2 (no sweep): t1, z0, z1, p1
    zf2 = build_family(("zf", 2, False))
    t1s = np.linspace(0.1, 0.9, 5 if args.quick else 17)
    zs = np.linspace(-3500, 3500, 7 if args.quick else 29)
    ps = np.arange(4 if args.quick else 8) / (4 if args.quick else 8)
    grid = [np.array([t, 1 - t, z0, z1, p]) for t in t1s for z0 in zs for z1 in zs for p in ps]
    chunks = np.array_split(np.arange(len(grid)), args.workers * 8)
    tasks += [(task_grid, (("zf", 2, False), [grid[i] for i in c])) for c in chunks]
    print(f"{len(tasks)} tasks ({len(grid)} grid points), {args.workers} workers", flush=True)

    # Optuna tasks are the longest: start them first.
    tasks.sort(key=lambda t: 0 if t[0] is task_optuna else 1)
    all_rows, results = [], []
    with Pool(args.workers, initializer=_init_worker) as pool:
        for i, res in enumerate(pool.imap_unordered(_dispatch, tasks)):
            all_rows += res.pop("rows")
            results.append(res)
            if (i + 1) % 20 == 0 or i + 1 == len(tasks):
                print(f"  {i + 1}/{len(tasks)} tasks, {len(all_rows)} evals, {time.time() - t0:.0f} s", flush=True)

        # local polishing of the best N = 2 grid points (with and without sweep)
        g = sorted([r for r in all_rows if r["source"] == "grid"], key=lambda r: -r["score"])
        seeds_x, seen = [], set()
        for r in g:
            key = tuple(np.round(json.loads(r["x"]), 3))
            if key not in seen:
                seen.add(key)
                seeds_x.append(np.array(json.loads(r["x"])))
            if len(seeds_x) >= (2 if args.quick else 10):
                break
        pol = []
        for k, x in enumerate(seeds_x):
            u = (x - zf2.lo) / (zf2.hi - zf2.lo)
            pol.append((task_cma, (("zf", 2, False), 1000 + k, int(1500 * bscale), u, 0.05, "grid+cma")))
            fs = build_family(("zf", 2, True))
            us = np.concatenate([u, np.full(2, 0.5)])  # sweeps start at 0
            pol.append((task_cma, (("zf", 2, True), 1000 + k, int(2000 * bscale), us, 0.05, "grid+cma")))
        for res in pool.imap_unordered(_dispatch, pol):
            all_rows += res.pop("rows")
            results.append(res)
    print(f"search done: {len(all_rows)} evaluations, {time.time() - t0:.0f} s", flush=True)

    # ---- log ---------------------------------------------------------------------
    with gzip.open(OUT / "evals.csv.gz", "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "family", "run", "x"] + LOG_KEYS)
        w.writeheader()
        w.writerows(all_rows)

    # ---- per configuration summary --------------------------------------------------
    configs = {}
    for res in results:
        if "best_x" not in res:
            continue
        name = build_family(res["spec"]).name
        c = configs.setdefault(name, {"spec": res["spec"], "runs": []})
        c["runs"].append((res["source"], res["seed"], res["best_score"], res["best_x"]))
    # the grid itself counts as a run for zoned_focus_2
    if g:
        configs["zoned_focus_2"]["runs"].append(("grid", 0, g[0]["score"], np.array(json.loads(g[0]["x"]))))

    ref = Evaluator(RadialPSF(n_pupil=2048, dr_nm=10.0), Targets(dz_nm=12.5))
    crlb_hi = LateralCRLB(Photometry(**{**HIGH_BACKGROUND.__dict__}))
    summary, curves = {}, {}
    for name, c in configs.items():
        fam = build_family(c["spec"])
        runs = sorted(c["runs"], key=lambda r: -r[2])
        cma_scores = sorted([r[2] for r in runs if r[0] == "cma"], reverse=True)
        best_src, _, best_search, bx = next((r for r in runs if r[0] != "optuna_tpe"), runs[0])
        m = ref.evaluate_mask(*fam.mask(bx), with_crlb=True, keep_curves=True)
        c_hi = crlb_hi.sigma_x(ref.rad.r_nm, m["curve_intensity"])
        entry = {
            "family": name, "dim": fam.dim, "best_source": best_src, "search_score": float(best_search),
            "x": [float(v) for v in bx], "mask": fam.describe(bx),
            "cma_restart_scores": [float(s) for s in cma_scores],
            "cma_hit_rate_1pct": float(np.mean(np.array(cma_scores) >= 0.99 * cma_scores[0])) if cma_scores else None,
            "crlb_max_highbg_nm": float(c_hi.max()), "crlb_mean_highbg_nm": float(c_hi.mean()),
        }
        entry.update({k: (float(v) if not isinstance(v, (bool, np.bool_)) else bool(v))
                      for k, v in m.items() if not k.startswith("curve")})
        opt = [r for r in runs if r[0] == "optuna_tpe"]
        if opt:
            mo = ref.evaluate_mask(*fam.mask(opt[0][3]))
            entry["optuna_score"] = float(mo["score"])
            entry["optuna_softmin"] = float(mo["ee_softmin"])
        summary[name] = entry
        m["curve_crlb_highbg_nm"] = c_hi
        m["phase_rho"] = ref.rad.rho
        m["phase_waves"] = ref.phase(*fam.mask(bx))
        curves[name] = m

    order = sorted(summary, key=lambda k: -summary[k]["score"])
    cols = ["score", "ee_softmin", "ee_min", "ee_mean", "eta", "ee_flatness", "strehl_min", "strehl_mean",
            "ee50_max_core_nm", "crlb_max_nm", "crlb_max_highbg_nm", "feasible"]
    print("\n## Best mask per configuration (reference evaluator, |z| <= 3 um)\n")
    print("| configuration | dim | " + " | ".join(cols) + " | CMA hit rate | Optuna |")
    print("|---|---|" + "---|" * (len(cols) + 2))
    for k in order:
        e = summary[k]
        vals = []
        for c in cols:
            v = e[c]
            vals.append(("yes" if v else "no") if isinstance(v, bool) else (f"{v:.0f}" if abs(v) >= 100 else f"{v:.3f}"))
        hr = "–" if e["cma_hit_rate_1pct"] is None else f"{e['cma_hit_rate_1pct']:.2f}"
        op = f"{e['optuna_score']:.3f}" if "optuna_score" in e else "–"
        print(f"| {k} | {e['dim']} | " + " | ".join(vals) + f" | {hr} | {op} |")
    print()
    for k in order:
        print(f"  {k}: {summary[k]['mask']}")

    js = {"targets": ref.describe(), "n_evals": len(all_rows), "configs": summary, "order": order,
          "runtime_s": time.time() - t0}
    (OUT / "results.json").write_text(json.dumps(js, indent=2, default=float))
    np.savez_compressed(OUT / "best_curves.npz", **{f"{k}__{kk}": v for k, m in curves.items()
                                                    for kk, v in m.items()
                                                    if (kk.startswith("curve") or kk.startswith("phase"))
                                                    and kk != "curve_intensity"})
    make_figures(summary, curves, order, ref)
    print(f"\ndone in {time.time() - t0:.0f} s; wrote {OUT}")


def _dispatch(t):
    fn, a = t
    return fn(a)


def make_figures(summary, curves, order, ref):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = OUT / "figures"
    ceil = ref.describe()["flat_ceiling"]
    base_line = {k: summary[k]["score"] for k in ("spherical+defocus", "zernike_m0") if k in summary}

    # 1. gate: best score vs number of zones
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for label, fmt, col in [("zoned_focus_{n}", "o-", "#2b6cb0"), ("zoned_focus_{n}_sweep", "s-", "#2f855a"),
                            ("layer_cake_{n}", "^-", "#805ad5"), ("annular_binary_{e}", "D-", "#d69e2e")]:
        ns = [n for n in (2, 3, 4, 5) if label.format(n=n, e=n - 1) in summary]
        es = [summary[label.format(n=n, e=n - 1)] for n in ns]
        feas = np.array([e["feasible"] for e in es])
        vals = np.array([e["ee_softmin"] for e in es])
        name = label.replace("_{n}", " N").replace("_{e}", " (N-1 edges)")
        ax.plot(np.array(ns)[feas], vals[feas], fmt, color=col, label=name)
        ax.plot(np.array(ns)[~feas], vals[~feas], fmt[0], color=col, mfc="none", ls="none")
        for n in ns:
            sc = [v for v in summary[label.format(n=n, e=n - 1)]["cma_restart_scores"] if v > 0]
            off = {"#2b6cb0": -0.06, "#2f855a": 0.06}.get(col, 0.0)
            if sc:
                ax.scatter([n + off] * len(sc), sc, s=8, color=col, alpha=0.3)
    for (k, v), ls in zip(base_line.items(), ("--", ":")):
        ax.axhline(v, color="#c05621", ls=ls, lw=1.2, label=f"{k} (re-optimised)")
    ax.axhline(ceil, color="k", ls="-.", lw=0.8, label=f"flat ceiling {ceil:.3f}")
    ax.axhspan(0.15, 0.19, color="#f6e05e", alpha=0.2, label="agreed target 0.15–0.19")
    ax.set_xticks([2, 3, 4, 5])
    ax.set_xlabel("number of zones N")
    ax.set_ylabel("soft-min EE400 over ±3 µm (S ≥ 0.03)")
    ax.set_ylim(0, ceil * 1.08)
    ax.set_title("Phase 2 gate: best mask vs N (small dots: CMA-ES restarts; hollow: infeasible, S < 0.03)",
                 fontsize=9)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "gate_vs_N.png", dpi=150)
    plt.close(fig)

    # 2. curves of the best masks
    best_zoned = max((k for k in summary if k.startswith("zoned_focus")), key=lambda k: summary[k]["score"])
    best_cake = max((k for k in summary if k.startswith("layer_cake")), key=lambda k: summary[k]["score"])
    show = [k for k in ("spherical+defocus", "zernike_m0", best_cake, "annular_binary_3", best_zoned)
            if k in curves]
    colors = ["#c05621", "#d69e2e", "#805ad5", "#718096", "#2b6cb0"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for col, k in zip(colors, show):
        m = curves[k]
        z = m["curve_z_nm"] / 1000
        lw = 2.2 if k == best_zoned else 1.4
        axes[0, 0].plot(z, m["curve_encircled_energy"], color=col, lw=lw, label=f"{k} ({summary[k]['score']:.3f})")
        axes[0, 1].plot(z, m["curve_strehl"], color=col, lw=lw)
        axes[1, 0].plot(z, m["curve_ee50_nm"], color=col, lw=lw)
        axes[1, 1].plot(z, m["curve_crlb_nm"], color=col, lw=lw)
        axes[1, 1].plot(z, m["curve_crlb_highbg_nm"], color=col, lw=0.9, ls="--")
    axes[0, 0].axhline(ceil, color="k", ls=":", lw=0.8)
    axes[0, 0].axhspan(0.15, 0.19, color="#f6e05e", alpha=0.2)
    axes[0, 1].axhline(0.03, color="k", ls=":", lw=0.8)
    labels = ["EE400", "Strehl (peak)", "EE50 (nm, reported only)",
              "lateral CRLB (nm): solid nominal, dashed high background"]
    for ax, lab in zip(axes.ravel(), labels):
        ax.set_ylabel(lab, fontsize=9)
        ax.grid(alpha=0.3)
    axes[1, 1].set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("z (µm)")
    axes[0, 0].legend(fontsize=7.5)
    fig.suptitle("Best mask of each family (reference evaluator); legend: soft-min EE400")
    fig.tight_layout()
    fig.savefig(fig_dir / "curves.png", dpi=150)
    plt.close(fig)

    # 3. phase profiles
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for col, k in zip(colors, show):
        ph = curves[k]["phase_waves"]
        ax.plot(curves[k]["phase_rho"], ph - ph[0], color=col, lw=2.0 if k == best_zoned else 1.2, label=k)
    ax.set_xlabel("normalised pupil radius ρ")
    ax.set_ylabel("mask phase (waves, unwrapped, centre = 0)")
    ax.set_title("Phase profiles of the best masks", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "phase_profiles.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
