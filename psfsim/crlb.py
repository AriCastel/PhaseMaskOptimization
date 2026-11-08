"""
Lateral Cramér-Rao lower bound for a camera-sampled, radially symmetric PSF.

Model per camera pixel k (Poisson signal + background, sCMOS read noise via
the usual Poisson-shift approximation):

    mu_k(x0)  = N * Integral_pixel I(|r - r0|) dA + b
    F_xx      = sum_k (d mu_k / d x0)^2 / (mu_k + sigma_read^2)
    sigma_x   = 1 / sqrt(F_xx)

I is the unit-power PSF from `RadialPSF` (1/nm^2), so N is the number of
detected photons in the *whole* PSF; photons falling outside the fitting ROI
are lost, as they would be in a real ROI-based fit.

For an emitter at a pixel centre the PSF is mirror-symmetric about the
emitter, so the Fisher cross terms between x and (y, z, N, b) vanish and
1/sqrt(F_xx) is the exact lateral bound even with z, N and b unknown. Off-centre
emitters can be included through `offsets_nm` (the bound is averaged in
variance over them, ignoring the then non-zero cross terms).

PLACEHOLDER photometry: the photon count, background and read noise defaults
below are generic single-molecule numbers, NOT ReflEx3D measurements. Replace
them with values from real frames before trusting absolute CRLB numbers (plan,
open question).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class Photometry:
    n_photons: float = 2000.0  # detected signal photons per frame (PLACEHOLDER)
    background: float = 10.0  # photons / pixel / frame (PLACEHOLDER)
    read_noise_e: float = 1.6  # sCMOS rms read noise, e- (PLACEHOLDER)
    camera_pixel_nm: float = 108.33  # object-side pitch (6.5 um / 60x)
    roi_px: int = 21  # fitting window, pixels per side (odd)
    subsample: int = 8  # sub-samples per pixel side for area integration


class LateralCRLB:
    """Precomputes the ROI sampling geometry; call with I(r, z) stacks."""

    def __init__(self, photometry: Photometry | None = None,
                 offsets_nm: Sequence[tuple[float, float]] = ((0.0, 0.0),)) -> None:
        self.ph = photometry or Photometry()
        p = self.ph
        if p.roi_px % 2 == 0:
            raise ValueError("roi_px must be odd (emitter on the central pixel)")
        s = p.subsample
        sub = ((np.arange(s) + 0.5) / s - 0.5) * p.camera_pixel_nm
        centres = (np.arange(p.roi_px) - p.roi_px // 2) * p.camera_pixel_nm
        coords = (centres[:, None] + sub[None, :]).ravel()  # (roi*s,)
        self._geom = []
        for ox, oy in offsets_nm:
            xx, yy = np.meshgrid(coords - ox, coords - oy, indexing="xy")
            r = np.hypot(xx, yy)
            cos_x = np.divide(xx, r, out=np.zeros_like(xx), where=r > 0)
            self._geom.append((r.ravel(), cos_x.ravel()))
        self._area = (p.camera_pixel_nm / s) ** 2
        self._shape = (p.roi_px, s, p.roi_px, s)

    def _pixel_sum(self, v: np.ndarray) -> np.ndarray:
        # v is on the (y-sub, x-sub) grid flattened as (roi*s, roi*s)
        return v.reshape(self._shape).sum(axis=(1, 3)) * self._area

    def sigma_x(self, r_nm: np.ndarray, intensity: np.ndarray) -> np.ndarray:
        """Lateral CRLB (nm) for each plane of `intensity` (n_z, n_r)."""
        p = self.ph
        inten = np.atleast_2d(intensity)
        dr = r_nm[1] - r_nm[0]
        dI = np.gradient(inten, dr, axis=1)
        out = np.empty(inten.shape[0])
        for k in range(inten.shape[0]):
            var = 0.0
            for r, cos_x in self._geom:
                i_s = np.interp(r, r_nm, inten[k], right=0.0)
                d_s = np.interp(r, r_nm, dI[k], right=0.0)
                mu = p.n_photons * self._pixel_sum(i_s) + p.background
                dmu = -p.n_photons * self._pixel_sum(d_s * cos_x)
                fisher = np.sum(dmu**2 / (mu + p.read_noise_e**2))
                var += 1.0 / fisher if fisher > 0 else np.inf
            out[k] = np.sqrt(var / len(self._geom))
        return out
