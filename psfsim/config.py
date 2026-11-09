"""
Optical system description shared by every stage of the pipeline.

NOTE (Phase 0): the original `config.py` was not part of the uploaded
skeleton. This is a reconstruction from how the other modules use
`OpticalSystem`. Defaults are the EDOF target system (NA 1.1, 525 nm,
water n = 1.33). If the lab copy differs (apodization, grid, SLM or camera
numbers), replace this file and re-run `scripts/phase0_check.py`.

Sampling
--------
The 2D model is a single FFT from the pupil to the image plane. With the pupil
radius spanning `pupil_radius_px` samples of a `grid_size` grid:

    sim_pixel_nm     = pupil_radius_px * wavelength / (grid_size * NA)
    field_of_view_nm = grid_size * sim_pixel_nm = pupil_radius_px * wavelength / NA

So a larger pupil radius (in samples) gives a larger field of view and a
coarser image grid; a larger grid at fixed pupil radius gives a finer image
grid. Defaults (1024, 128) give ~60 nm sampling over a ~61 um field.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

APODIZATIONS = ("emission", "uniform", "illumination")


@dataclass(frozen=True)
class OpticalSystem:
    # Objective / sample
    na: float = 1.1
    wavelength_nm: float = 525.0
    n_medium: float = 1.33

    # Pupil amplitude model (see pupil.pupil_amplitude):
    #   "uniform"      flat pupil (default: reproduces the existing
    #                  tracking-metrics figure, e.g. EE400(0) = 0.88 no mask)
    #   "emission"     1/sqrt(cos theta): isotropic emitter collected by an
    #                  aplanatic objective (physical detection PSF; lowers the
    #                  EE budget by ~10 %, see docs/phase0_report.md)
    #   "illumination" sqrt(cos theta): aplanatic focusing of a flat beam
    apodization: str = "uniform"

    # Simulation grid
    grid_size: int = 1024
    pupil_radius_px: float = 128.0

    # Camera (object-side pitch = camera_pixel_um / magnification)
    camera_pixel_um: float = 6.5
    magnification: float = 66.0  # 6.5 um / 66 = 98.5 nm on the sample

    # SLM, as imaged onto the pupil
    slm_pixels_across_pupil: float = 600.0  # SLM pixels spanning the pupil diameter
    slm_bit_depth: Optional[int] = 8  # None = continuous phase levels
    slm_wrap: bool = True  # phase displayed modulo one wave

    def __post_init__(self) -> None:
        if not 0 < self.na < self.n_medium:
            raise ValueError("need 0 < NA < n_medium")
        if self.apodization not in APODIZATIONS:
            raise ValueError(f"apodization must be one of {APODIZATIONS}")

    # -- derived quantities -------------------------------------------------
    @property
    def sim_pixel_nm(self) -> float:
        return self.pupil_radius_px * self.wavelength_nm / (self.grid_size * self.na)

    @property
    def field_of_view_nm(self) -> float:
        return self.grid_size * self.sim_pixel_nm

    @property
    def camera_pixel_object_nm(self) -> float:
        return self.camera_pixel_um * 1000.0 / self.magnification

    @property
    def cutoff_frequency(self) -> float:
        """Incoherent cutoff 2 NA / lambda, in cycles per nm."""
        return 2.0 * self.na / self.wavelength_nm

    @property
    def depth_of_field_nm(self) -> float:
        """Native depth of field n * lambda / NA^2."""
        return self.n_medium * self.wavelength_nm / self.na**2

    @property
    def sin_max(self) -> float:
        return self.na / self.n_medium

    @property
    def cos_min(self) -> float:
        """c0 = cos(theta_max), the lower limit of the pupil in c = cos(theta)."""
        return float(np.sqrt(1.0 - self.sin_max**2))

    def with_(self, **changes) -> "OpticalSystem":
        """Copy with some fields changed."""
        return replace(self, **changes)
