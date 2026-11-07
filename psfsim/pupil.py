"""
Pupil grid, amplitude and depth-induced defocus.

NOTE (Phase 0): reconstructed; the original module was not in the skeleton.
The functions here only use `grid.rho` and `grid.mask`, so they work both on
the 2D FFT grid and on the 1D radial grid used by `psfsim.radial`. That is
what guarantees the two models use the same apodization and defocus.

Conventions
-----------
* rho is the normalised pupil radius, rho = sin(theta) / sin(theta_max).
* Phases are in waves.
* Defocus: a source at depth z adds  (z / lambda) * (n cos(theta) - n)  waves
  (exact), or its paraxial limit  -(z / lambda) * NA^2 rho^2 / (2 n). The
  piston n z / lambda is dropped; it does not change the intensity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import OpticalSystem


@dataclass
class PupilGrid:
    x: np.ndarray
    y: np.ndarray
    rho: np.ndarray
    theta: np.ndarray
    mask: np.ndarray  # True inside the pupil
    pupil_radius_px: float = float("nan")

    @property
    def n_inside(self) -> int:
        return int(self.mask.sum())


def make_pupil_grid(system: OpticalSystem) -> PupilGrid:
    """Square grid whose zero frequency sits at index grid_size // 2."""
    n = system.grid_size
    r = system.pupil_radius_px
    idx = (np.arange(n) - n // 2) / r
    x, y = np.meshgrid(idx, idx, indexing="xy")
    rho = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return PupilGrid(x=x, y=y, rho=rho, theta=theta, mask=rho <= 1.0, pupil_radius_px=r)


def radial_pupil_grid(rho: np.ndarray) -> PupilGrid:
    """A 1D 'grid' of radii along theta = 0, for radially symmetric evaluation."""
    rho = np.asarray(rho, dtype=float)
    return PupilGrid(
        x=rho.copy(),
        y=np.zeros_like(rho),
        rho=rho,
        theta=np.zeros_like(rho),
        mask=rho <= 1.0,
    )


def cos_theta(grid: PupilGrid, system: OpticalSystem) -> np.ndarray:
    s = np.clip(system.sin_max * grid.rho, 0.0, 1.0)
    return np.sqrt(1.0 - s**2)


def pupil_amplitude(grid: PupilGrid, system: OpticalSystem) -> np.ndarray:
    c = np.where(grid.mask, cos_theta(grid, system), 1.0)
    if system.apodization == "emission":
        a = 1.0 / np.sqrt(c)
    elif system.apodization == "illumination":
        a = np.sqrt(c)
    else:
        a = np.ones_like(c)
    return np.where(grid.mask, a, 0.0)


def defocus_phase(
    grid: PupilGrid, system: OpticalSystem, z_nm: float, exact: bool = True
) -> np.ndarray:
    """Depth-induced phase in waves for a source at z_nm (positive z = the
    sign convention used everywhere in psfsim)."""
    n = system.n_medium
    if exact:
        w = (z_nm / system.wavelength_nm) * n * (cos_theta(grid, system) - 1.0)
    else:
        w = -(z_nm / system.wavelength_nm) * system.na**2 * grid.rho**2 / (2.0 * n)
    return np.where(grid.mask, w, 0.0)


def _inside(phase: np.ndarray, grid: PupilGrid) -> np.ndarray:
    return np.asarray(phase)[grid.mask]


def rms_over_pupil(phase: np.ndarray, grid: PupilGrid) -> float:
    v = _inside(phase, grid)
    return float(np.sqrt(np.mean((v - v.mean()) ** 2))) if v.size else 0.0


def peak_to_valley(phase: np.ndarray, grid: PupilGrid) -> float:
    v = _inside(phase, grid)
    return float(v.max() - v.min()) if v.size else 0.0
