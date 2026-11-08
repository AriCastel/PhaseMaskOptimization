"""
Phase 1: the scalar objective and constraints for EDOF mask design.

Primary objective (plan, Phase 1)
---------------------------------
    maximise   softmin_{|z| <= 3 um} EE400(z)
    subject to S(z)   >= S_min   for all |z| <= 3 um     (peak Strehl, detection proxy)
               EE50(z) <= R_max  for all |z| <= 2 um     (spot-size cap near focus)

Reported alongside: eta (in-range share of the mask-independent P(r_opt)
budget), flatness (min/mean), CV, mean Strehl, mean PW-M2 and, optionally,
the worst-case lateral CRLB.

Soft minimum
------------
A power mean with negative exponent p,  M_p(v) = (mean v^p)^(1/p),  p < 0.
It is scale-free (no temperature in EE units), lies between min(v) (p -> -inf)
and the harmonic mean (p = -1), and is smooth, so gradient and CMA-type
optimisers see progress before the single worst plane moves. p = -10 is the
default: a plane at half the level of the rest is weighted ~1000x more.

Thresholds
----------
S_min and R_max are PROVISIONAL until agreed with the advisor (Phase 1 gate).
`scripts/phase1_pareto.py` shows how the attainable EE floor depends on them.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping

import numpy as np

from .crlb import LateralCRLB, Photometry
from .radial import RadialPSF


def soft_min(values: np.ndarray, p: float = -10.0) -> float:
    v = np.asarray(values, dtype=float)
    if np.any(v <= 0):
        return float(v.min())
    return float(np.mean(v**p) ** (1.0 / p))


@dataclass(frozen=True)
class Targets:
    z_range_nm: tuple[float, float] = (-3000.0, 3000.0)
    dz_nm: float = 50.0
    # PROVISIONAL until agreed with the advisor (Phase 1 gate). Chosen at the
    # knee of the threshold map in docs/phase1_report.md: tightening either
    # one further costs EE floor; 1500 nm ~ the no-mask EE50 at |z| = 2 um.
    s_min: float = 0.03
    r_max_nm: float = 1500.0
    ee50_range_nm: tuple[float, float] = (-2000.0, 2000.0)
    softmin_p: float = -10.0
    penalty: float = 2.0  # weight on relative constraint violation in `score`

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
            "ee_softmin": soft_min(ee, t.softmin_p),
            "ee_min": float(ee.min()),
            "ee_mean": float(ee.mean()),
            "ee_flatness": float(ee.min() / ee.mean()),
            "ee_cv": float(ee.std() / ee.mean()),
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
        out["score"] = out["ee_softmin"] - t.penalty * out["violation"]
        if with_crlb:
            c = self.crlb.sigma_x(self.rad.r_nm, inten)
            out["crlb_max_nm"] = float(c.max())
            out["crlb_mean_nm"] = float(c.mean())
            if keep_curves:
                out["curve_crlb_nm"] = c
        if keep_curves:
            out.update({f"curve_{k}": v for k, v in m.items()})
            out["curve_z_nm"] = self.z
        return out

    def evaluate_mask(self, mask=None, mask_params=None, **kw) -> dict:
        return self.evaluate(self.phase(mask, mask_params), **kw)

    def describe(self) -> dict:
        return {"targets": asdict(self.targets), "photometry": asdict(self.crlb.ph),
                "budget_power_ropt_nm": self.budget_p,
                "flat_ceiling": self.budget_p / (self.targets.z_range_nm[1] - self.targets.z_range_nm[0])}
