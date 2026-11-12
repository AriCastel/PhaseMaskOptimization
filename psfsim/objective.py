"""
Phase 1: the scalar objective and constraints for EDOF mask design.

Primary objective (plan, Phase 1)
---------------------------------
    maximise   softmin_{|z| <= 3 um} EE400(z)
    subject to S(z)   >= S_min   for all |z| <= 3 um     (peak Strehl, detection proxy)

EE50 is NOT a constraint (decision 2026-11-08): beyond ~2 um it tracks the
geometric defocus cone and, since EE400 <= ~0.28, it mostly measures the halo.
It is reported (`ee50_max_core_nm`, max over |z| <= 2 um) for a-posteriori
interpretation. An R_max can still be set for experiments.

Reported alongside: eta (in-range share of the mask-independent P(r_opt)
budget), flatness (min/mean), CV, mean Strehl, mean PW-M2 and, optionally,
the worst-case lateral CRLB.

Soft minimum
------------
A power mean with negative exponent p,  M_p(v) = (mean_z v^p)^(1/p),  p < 0,
with the mean taken as a trapezoid-weighted z-integral (grid independent).
It is scale-free (no temperature in EE units), lies between min(v) (p -> -inf)
and the harmonic mean (p = -1), and is smooth, so gradient and CMA-type
optimisers see progress before the single worst plane moves. p = -10 is the
default: a plane at half the level of the rest is weighted ~1000x more.

EE denominator
--------------
EE400 = P(r_opt) / P(total_radius = 10 um) is a ratio inside a finite window, so
a mask can raise it by pushing light OUT of the window (a zone defocused by more
than ~7 um spreads its light beyond 10 um): the denominator shrinks, the core
does not gain. The Phase 2 families could not do this (zone depths bounded at
+-4 um); families with free per-zone Zernikes can, and CMA-ES finds it (see
docs/zoned_zernike6_report.md). `Targets.ee_denominator = "emitted"` scores
P(r_opt) as a fraction of ALL the light through the pupil instead, which cannot
be gamed; "window" (default) is the agreed Phase 1-2 definition. Both soft-mins
are always reported (`ee_softmin`: window, `p_softmin`: emitted), together with
`window_fraction_min` = min_z P(10 um), which is ~0.95-0.99 for a mask that keeps
its light.

Thresholds
----------
S_min = 0.03, no EE50 cap, target EE400 floor 0.15-0.19 (agreed 2026-11-08;
see docs/phase1_report.md for the threshold map behind them).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping

import numpy as np

from .crlb import LateralCRLB, Photometry
from .radial import RadialPSF


def soft_min(values: np.ndarray, p: float = -10.0, weights: np.ndarray | None = None) -> float:
    """Power mean (sum w v^p / sum w)^(1/p). Pass trapezoid weights in z so the
    value is a z-integral mean and does not depend on the grid spacing (with
    p = -10 the range ends, where EE usually drops, otherwise get a
    grid-dependent weight)."""
    v = np.asarray(values, dtype=float)
    if np.any(v <= 0):
        return float(v.min())
    w = np.ones_like(v) if weights is None else np.asarray(weights, dtype=float)
    return float((np.sum(w * v**p) / np.sum(w)) ** (1.0 / p))


def trapezoid_weights(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    w = np.zeros_like(z)
    d = np.diff(z)
    w[:-1] += d / 2
    w[1:] += d / 2
    return w


@dataclass(frozen=True)
class Targets:
    z_range_nm: tuple[float, float] = (-3000.0, 3000.0)
    dz_nm: float = 50.0
    # Agreed 2026-11-08 (Phase 1 gate): Strehl floor at the knee of the
    # threshold map; EE50 reported only (r_max_nm = inf); EE400 target 0.15-0.19.
    s_min: float = 0.03
    r_max_nm: float = float("inf")
    ee50_range_nm: tuple[float, float] = (-2000.0, 2000.0)
    softmin_p: float = -10.0
    penalty: float = 2.0  # weight on relative constraint violation in `score`
    # What the scored encircled energy is a fraction of: "window" = the light inside
    # total_radius (EE400 as agreed), "emitted" = all the light through the pupil.
    ee_denominator: str = "window"

    def __post_init__(self) -> None:
        if self.ee_denominator not in ("window", "emitted"):
            raise ValueError("ee_denominator must be 'window' or 'emitted'")

    @property
    def z_nm(self) -> np.ndarray:
        lo, hi = self.z_range_nm
        n = int(round((hi - lo) / self.dz_nm))
        return np.linspace(lo, hi, n + 1)


class Evaluator:
    """Scores a radially symmetric pupil phase against `Targets`."""

    def __init__(self, rad: RadialPSF | None = None, targets: Targets | None = None,
                 photometry: Photometry | None = None) -> None:
        self.rad = rad or RadialPSF()
        self.targets = targets or Targets()
        self.z = self.targets.z_nm
        self._w = trapezoid_weights(self.z)
        self.budget_p = self.rad.budgets().power_ropt_nm
        self.crlb = LateralCRLB(photometry or Photometry(camera_pixel_nm=self.rad.system.camera_pixel_object_nm))
        lo, hi = self.targets.ee50_range_nm
        self._core = (self.z >= lo - 1e-9) & (self.z <= hi + 1e-9)
        self.n_evals = 0

    # -- phase helpers ------------------------------------------------------
    def phase(self, mask=None, mask_params: Mapping | None = None, aberrations=None,
              phase_waves=None) -> np.ndarray:
        return self.rad.total_phase(mask, mask_params, aberrations, phase_waves)

    # -- evaluation ---------------------------------------------------------
    def evaluate(self, phase: np.ndarray, with_crlb: bool = False,
                 keep_curves: bool = False) -> dict:
        t = self.targets
        m, inten = self.rad._metrics(self.z, phase, self.rad.reference_power)
        self.n_evals += 1
        ee, s, p = m["encircled_energy"], m["strehl"], m["power_ropt"]
        ee50_core = m["ee50_nm"][self._core]
        out = {
            "ee_softmin": soft_min(ee, t.softmin_p, self._w),
            "ee_min": float(ee.min()),
            "ee_mean": float(ee.mean()),
            "ee_flatness": float(ee.min() / ee.mean()),
            "ee_cv": float(ee.std() / ee.mean()),
            "p_softmin": soft_min(p, t.softmin_p, self._w),
            "p_min": float(p.min()),
            "p_mean": float(p.mean()),
            "window_fraction_min": float(m["total_power"].min()),
            "eta": float(np.trapezoid(p, self.z) / self.budget_p),
            "strehl_min": float(s.min()),
            "strehl_mean": float(s.mean()),
            "ee50_max_core_nm": float(np.nanmax(ee50_core)) if np.isfinite(ee50_core).any() else np.inf,
            "pwm2_mean_nm2": float(np.mean(m["power_weighted_m2"])),
        }
        v_s = max(0.0, (t.s_min - out["strehl_min"]) / t.s_min) if t.s_min > 0 else 0.0
        v_r = (max(0.0, (out["ee50_max_core_nm"] - t.r_max_nm) / t.r_max_nm)
               if np.isfinite(t.r_max_nm) else 0.0)
        out["violation"] = v_s + v_r
        out["feasible"] = out["violation"] == 0.0
        objective = out["p_softmin"] if t.ee_denominator == "emitted" else out["ee_softmin"]
        out["score"] = objective - t.penalty * out["violation"]
        if with_crlb:
            c = self.crlb.sigma_x(self.rad.r_nm, inten)
            out["crlb_max_nm"] = float(c.max())
            out["crlb_mean_nm"] = float(c.mean())
            if keep_curves:
                out["curve_crlb_nm"] = c
        if keep_curves:
            out.update({f"curve_{k}": v for k, v in m.items()})
            out["curve_intensity"] = inten
            out["curve_z_nm"] = self.z
        return out

    def evaluate_mask(self, mask=None, mask_params=None, **kw) -> dict:
        return self.evaluate(self.phase(mask, mask_params), **kw)

    def describe(self) -> dict:
        return {"targets": asdict(self.targets), "photometry": asdict(self.crlb.ph),
                "budget_power_ropt_nm": self.budget_p,
                "flat_ceiling": self.budget_p / (self.targets.z_range_nm[1] - self.targets.z_range_nm[0])}
