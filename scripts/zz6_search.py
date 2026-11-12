"""
Zoned radial-Zernike masks: free Z(n, 0), n <= 6, in every zone (N = 3 and 4),
against the `zoned_focus` (exact defocus + piston [+ sweep]) candidates of Phase 2.

Objective: maximise the soft-min over |z| <= 3 um of the energy within r = 400 nm,
subject to peak Strehl >= 0.03. `--denominator window` takes that energy as a
fraction of the light inside 10 um (EE400 exactly as in Phase 2); `--denominator
emitted` (default) as a fraction of all the light through the pupil, which a mask
cannot inflate by throwing light out of the 10 um window (see psfsim/objective.py).
Search evaluator: 1024 pupil
rings, dr 20 nm, dz 50 nm. Final scores: reference evaluator (2048 rings, dr
10 nm, dz 12.5 nm) with the lateral CRLB, nominal and high background.

Search stages for the new family `zoned_zernike6_N` (psfsim.families.ZonedZernikeRadial):
  cold        20 random CMA-ES restarts, the Phase 2 protocol (sigma 0.3, budget 250 x dim <= 4000)
  cold-long   12 random restarts, 3x population, 12 000 evaluations
  warm        every distinct Phase 2 run of zoned_focus_N and zoned_focus_N + sweep, projected
              onto the family (least squares per zone), CMA-ES sigma 0.05 and 0.15, 6000 evaluations
  polish      the 6 best distinct runs of the above, CMA-ES sigma 0.03, 15 000 evaluations
Matched-effort control for zoned_focus_N and zoned_focus_N + sweep (N = 3, 4), so that
the comparison is not new-family-with-more-search against old-family-with-less:
  cold-long   8 random restarts, 3x population, 12 000 evaluations
  polish      every distinct Phase 2 run, CMA-ES sigma 0.05, 10 000 evaluations
Order ablation: the best order-6 mask truncated to order 4 and to order 2
(paraxial defocus + piston) and re-polished (sigma 0.05, 8000 evaluations, 2 seeds).

Writes docs/zoned_zernike6/<denominator>/: results.json, best_curves.npz,
evals_top_per_family.csv.gz (the 2000 best evaluations per family).

Usage:  python scripts/zz6_search.py [--workers 4] [--quick]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import gzip  # noqa: E402
import heapq  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import RadialPSF  # noqa: E402
from psfsim.crlb import HIGH_BACKGROUND, LateralCRLB, Photometry  # noqa: E402
from psfsim.families import SphericalDefocus, ZonedFocus, ZonedZernikeRadial  # noqa: E402
from psfsim.objective import Evaluator, Targets  # noqa: E402

BASE = ROOT / "docs" / "zoned_zernike6"
OUT = BASE / "window"  # set in main()
DENOM = "window"
PHASE2 = ROOT / "docs" / "phase2"
LOG_KEYS = ["score", "ee_softmin", "ee_min", "p_softmin", "p_min", "window_fraction_min", "eta", "strehl_min",
            "strehl_mean", "ee50_max_core_nm", "feasible"]
KEEP_ROWS = 300  # best evaluations kept per run for the log

_EV = None


def _init_worker(denom):
    global _EV
    _EV = Evaluator(RadialPSF(n_pupil=1024, dr_nm=20.0), Targets(ee_denominator=denom))


def build_family(spec):
    kind, *args = spec
    return {"sph": SphericalDefocus, "zf": ZonedFocus, "zz": ZonedZernikeRadial}[kind](*args)


def task_cma(args):
    import cma

    spec, seed, budget, x0_u, sigma0, source, pop_factor = args
    fam = build_family(spec)
    rng = np.random.default_rng(seed)
    u0 = rng.uniform(0.05, 0.95, fam.dim) if x0_u is None else np.clip(np.asarray(x0_u), 0, 1)
    span = fam.hi - fam.lo
    opts = {"bounds": [0, 1], "seed": int(seed) + 1, "maxfevals": budget, "verbose": -9,
            "tolfun": 1e-7, "tolx": 1e-7}
    if pop_factor != 1:
        opts["popsize"] = int(pop_factor * (4 + int(3 * np.log(fam.dim))))
    es = cma.CMAEvolutionStrategy(u0, sigma0, opts)
    heap, best, n = [], (-np.inf, None), 0
    if x0_u is not None:  # a warm start is itself a candidate
        x = fam.lo + u0 * span
        best, n = (_EV.evaluate_mask(*fam.mask(x))["score"], x), 1
    while not es.stop():
        U = es.ask()
        F = []
        for u in U:
            x = fam.lo + np.clip(u, 0, 1) * span
            m = _EV.evaluate_mask(*fam.mask(x))
            n += 1
            F.append(-m["score"])
            item = (m["score"], n, [round(float(v), 6) for v in x],
                    [bool(m[k]) if k == "feasible" else round(float(m[k]), 6) for k in LOG_KEYS])
            if len(heap) < KEEP_ROWS:
                heapq.heappush(heap, item)
            elif item[0] > heap[0][0]:
                heapq.heapreplace(heap, item)
            if m["score"] > best[0]:
                best = (m["score"], x)
        es.tell(U, F)
    rows = [dict(source=source, family=fam.name, run=seed, x=json.dumps(x), **dict(zip(LOG_KEYS, vals)))
            for _, _, x, vals in heap]
    return {"spec": spec, "family": fam.name, "source": source, "seed": seed, "best_score": float(best[0]),
            "best_x": [float(v) for v in best[1]], "n_evals": n, "rows": rows}


def budget_for(dim):  # Phase 2 rule
    return int(np.clip(250 * dim, 1000, 4000))


def phase2_runs(family_name):
    """Best point of every distinct run of a family in the Phase 2 log (best first)."""
    best = {}
    with gzip.open(PHASE2 / "evals_top2000_per_family.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if r["family"] != family_name:
                continue
            k, s = (r["source"], r["run"]), float(r["score"])
            if k not in best or s > best[k][0]:
                best[k] = (s, np.array(json.loads(r["x"])))
    return [v for v in sorted(best.values(), key=lambda v: -v[0])]


def distinct_best(results, k):
    """The k best runs whose optima differ (score to 1e-4 or parameters)."""
    out = []
    for r in sorted(results, key=lambda r: -r["best_score"]):
        if all(abs(r["best_score"] - o["best_score"]) > 1e-4 or
               not np.allclose(r["best_x"], o["best_x"], atol=5e-3) for o in out):
            out.append(r)
        if len(out) == k:
            break
    return out


def run_pool(pool, tasks, sink, label, t0):
    tasks = sorted(tasks, key=lambda a: -a[2])  # longest first
    for i, res in enumerate(pool.imap_unordered(task_cma, tasks)):
        sink.append(res)
        if (i + 1) % 10 == 0 or i + 1 == len(tasks):
            print(f"  [{label}] {i + 1}/{len(tasks)} runs, {sum(r['n_evals'] for r in sink)} evals, "
                  f"{time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--quick", action="store_true", help="tiny run for testing the pipeline")
    ap.add_argument("--denominator", choices=("window", "emitted"), default="emitted",
                    help="EE400 as a fraction of the light inside 10 um (agreed definition) or of all emitted light")
    args = ap.parse_args()
    global OUT, DENOM
    DENOM = args.denominator
    OUT = BASE / (DENOM + ("_quick" if args.quick else ""))
    OUT.mkdir(parents=True, exist_ok=True)
    q = args.quick
    b = 0.03 if q else 1.0
    t0 = time.time()
    search = Evaluator(RadialPSF(n_pupil=1024, dr_nm=20.0), Targets(ee_denominator=DENOM))

    new_specs = [("zz", 3, 6), ("zz", 4, 6)]
    ctl_specs = [("zf", n, sw) for n in (3, 4) for sw in (False, True)]
    tasks = []
    for spec in new_specs:
        fam = build_family(spec)
        n = spec[1]
        tasks += [(spec, s, int(budget_for(fam.dim) * b), None, 0.3, "cold", 1) for s in range(2 if q else 20)]
        tasks += [(spec, 100 + s, int(12000 * b), None, 0.3, "cold-long", 3) for s in range(2 if q else 12)]
        k = 0
        for sw in (False, True):
            src = ZonedFocus(n, sw)
            for _, x in phase2_runs(src.name)[: 1 if q else None]:
                m, p = src.mask(x)
                u = (fam.from_phase(p["radii"], search.rad.rho, search.phase(m, p)) - fam.lo) / (fam.hi - fam.lo)
                for sig in (0.05, 0.15):
                    tasks.append((spec, 200 + k, int(6000 * b), list(u), sig, f"warm:{src.name}", 1))
                    k += 1
        # the optima of the other objective are starting points too
        other = BASE / ("window" if DENOM == "emitted" else "emitted") / "results.json"
        if other.exists() and not q:
            runs = json.loads(other.read_text())["configs"][fam.name]["runs"]
            for r in distinct_best([{"best_score": r["score"], "best_x": r["x"]} for r in runs], 3):
                u = (np.array(r["best_x"]) - fam.lo) / (fam.hi - fam.lo)
                tasks.append((spec, 200 + k, int(6000 * b), list(u), 0.1, "warm:other-objective", 1))
                k += 1
    for spec in ctl_specs:
        fam = build_family(spec)
        tasks += [(spec, 100 + s, int(12000 * b), None, 0.3, "cold-long", 3) for s in range(2 if q else 8)]
        for k, (_, x) in enumerate(phase2_runs(fam.name)[: 1 if q else None]):
            u = (x - fam.lo) / (fam.hi - fam.lo)
            tasks.append((spec, 300 + k, int(10000 * b), list(u), 0.05, "polish", 1))
    print(f"stage 1: {len(tasks)} CMA-ES runs, {args.workers} workers", flush=True)

    results = []
    with Pool(args.workers, initializer=_init_worker, initargs=(DENOM,)) as pool:
        run_pool(pool, tasks, results, "stage 1", t0)

        # polish the best distinct runs of the new family
        pol = []
        for spec in new_specs:
            fam = build_family(spec)
            mine = [r for r in results if r["family"] == fam.name]
            for k, r in enumerate(distinct_best(mine, 2 if q else 6)):
                u = (np.array(r["best_x"]) - fam.lo) / (fam.hi - fam.lo)
                pol.append((spec, 400 + k, int(15000 * b), list(u), 0.03, "polish", 1))
        run_pool(pool, pol, results, "polish", t0)

        # order ablation: truncate the best order-6 mask and re-polish
        abl = []
        for spec in new_specs:
            fam = build_family(spec)
            bx = np.array(max((r for r in results if r["family"] == fam.name),
                              key=lambda r: r["best_score"])["best_x"])
            for order in (2, 4):
                sub = ZonedZernikeRadial(spec[1], order)
                x = np.array([bx[fam.names.index(nm)] for nm in sub.names])
                u = (x - sub.lo) / (sub.hi - sub.lo)
                abl += [(("zz", spec[1], order), 500 + s, int(8000 * b), list(u), 0.05, "truncated+polish", 1)
                        for s in range(2)]
        run_pool(pool, abl, results, "ablation", t0)
    n_evals = sum(r["n_evals"] for r in results)
    print(f"search done: {n_evals} evaluations, {time.time() - t0:.0f} s", flush=True)

    # ---- log: best evaluations per family --------------------------------------
    rows = {}
    for r in results:
        rows.setdefault(r["family"], []).extend(r.pop("rows"))
    with gzip.open(OUT / "evals_top_per_family.csv.gz", "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "family", "run", "x"] + LOG_KEYS)
        w.writeheader()
        for fam_rows in rows.values():
            w.writerows(sorted(fam_rows, key=lambda r: -r["score"])[:2000])

    # ---- reference re-scoring ------------------------------------------------------
    ref = Evaluator(RadialPSF(n_pupil=2048, dr_nm=10.0), Targets(dz_nm=12.5, ee_denominator=DENOM))
    crlb_hi = LateralCRLB(Photometry(**{**HIGH_BACKGROUND.__dict__}))
    p2 = json.loads((PHASE2 / "results.json").read_text())["configs"]

    def full(fam, x):
        m = ref.evaluate_mask(*fam.mask(x), with_crlb=True, keep_curves=True)
        c_hi = crlb_hi.sigma_x(ref.rad.r_nm, m["curve_intensity"])
        e = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in m.items()
             if not k.startswith("curve")}
        e.update({"crlb_max_highbg_nm": float(c_hi.max()), "crlb_mean_highbg_nm": float(c_hi.mean()),
                  "family": fam.name, "dim": fam.dim, "x": [float(v) for v in x], "mask": fam.describe(x)})
        cur = {k: v for k, v in m.items() if k.startswith("curve") and k != "curve_intensity"}
        cur.update({"curve_crlb_highbg_nm": c_hi, "phase_rho": ref.rad.rho, "phase_waves": ref.phase(*fam.mask(x))})
        return e, cur

    summary, curves = {}, {}
    by_fam = {}
    for r in results:
        by_fam.setdefault(r["family"], []).append(r)
    for name, runs in by_fam.items():
        fam = build_family(runs[0]["spec"])
        # the search evaluator is coarser than the reference: re-score the best distinct
        # runs and keep the best one that is still feasible on the reference
        cands = [full(fam, np.array(r["best_x"])) + (r,) for r in distinct_best(runs, 8)]
        if name in p2:  # the Phase 2 optimum stays a candidate for its own family
            e2, c2 = full(fam, np.array(p2[name]["x"]))
            cands.append((e2, c2, {"source": "phase2", "best_score": float(
                search.evaluate_mask(*fam.mask(p2[name]["x"]))["score"])}))
        e, cur, r = max(cands, key=lambda c: c[0]["score"])
        e.update({"best_source": r["source"], "search_score": r["best_score"],
                  "runs": [{"source": x["source"], "seed": x["seed"], "score": x["best_score"],
                            "n_evals": x["n_evals"], "x": x["best_x"]} for x in
                           sorted(runs, key=lambda x: -x["best_score"])]})
        if name in p2:  # the Phase 2 optimum under this run's objective
            e["phase2_score"] = e2["score"]
        summary[name], curves[name] = e, cur
    # references re-scored here so that every curve comes from the same evaluator
    for name, fam, x in [("spherical+defocus", SphericalDefocus(), p2["spherical+defocus"]["x"]),
                         ("zoned_focus_4 (Phase 2)", ZonedFocus(4), p2["zoned_focus_4"]["x"]),
                         ("zoned_focus_4_sweep (Phase 2)", ZonedFocus(4, True), p2["zoned_focus_4_sweep"]["x"])]:
        summary[name], curves[name] = full(fam, np.array(x))
    m = ref.evaluate(ref.phase(None), with_crlb=True, keep_curves=True)
    c_hi = crlb_hi.sigma_x(ref.rad.r_nm, m["curve_intensity"])
    summary["no mask"] = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in m.items()
                          if not k.startswith("curve")}
    summary["no mask"].update({"crlb_max_highbg_nm": float(c_hi.max()), "family": "none", "dim": 0})
    curves["no mask"] = {k: v for k, v in m.items() if k.startswith("curve") and k != "curve_intensity"}
    curves["no mask"].update({"curve_crlb_highbg_nm": c_hi, "phase_rho": ref.rad.rho,
                              "phase_waves": np.zeros(ref.rad.n_pupil)})

    cols = ["score", "p_softmin", "p_min", "ee_softmin", "ee_min", "window_fraction_min", "eta", "strehl_min",
            "ee50_max_core_nm", "crlb_max_nm", "crlb_max_highbg_nm", "feasible"]
    order = sorted(summary, key=lambda k: -summary[k]["score"])
    print("\n| configuration | dim | " + " | ".join(cols) + " | found by |")
    print("|---|---|" + "---|" * (len(cols) + 1))
    for k in order:
        e = summary[k]
        vals = [("yes" if e[c] else "no") if isinstance(e[c], bool) else
                (f"{e[c]:.0f}" if abs(e[c]) >= 100 else f"{e[c]:.4f}") for c in cols]
        print(f"| {k} | {e['dim']} | " + " | ".join(vals) + f" | {e.get('best_source', '–')} |")
    print()
    for k in order:
        if "mask" in summary[k]:
            print(f"  {k}: {summary[k]['mask']}")

    js = {"targets": ref.describe(), "ee_denominator": DENOM, "n_evals": n_evals, "configs": summary, "order": order,
          "runtime_s": time.time() - t0}
    (OUT / "results.json").write_text(json.dumps(js, indent=2, default=float))
    np.savez_compressed(OUT / "best_curves.npz", **{f"{k}__{kk}": v for k, c in curves.items() for kk, v in c.items()})
    print(f"\ndone in {time.time() - t0:.0f} s; wrote {OUT}")


if __name__ == "__main__":
    main()
