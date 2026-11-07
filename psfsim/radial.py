"""
Fast radial (1D Debye) model for radially symmetric pupils.

For a pupil phase that depends only on rho, the 2D Fourier transform of
`PSFSimulator` reduces to a Hankel transform. With u = rho^2 (equal-area
sampling of the pupil, so each node carries the same pupil area):

    U(r, z) = 2 pi (NA/lambda)^2 * Integral_0^1  A(rho) exp(2 pi i psi(rho, z))
                                       J0(2 pi NA rho r / lambda)  du / 2
    psi     = mask(rho) + aberration(rho) + z * D(rho)
    D(rho)  = (n / lambda) (cos theta - 1)      exact defocus, waves per nm

Changing variable to c = cos(theta) turns the z-dependence into a pure Fourier
kernel exp(i k n z c), which is what makes the depth budgets below exact.
A(rho) and D(rho) come from `psfsim.pupil` (the same functions the 2D model
uses), so apodization and defocus are identical in both models.

Normalisation: I(r, z) = |U|^2 / (pupil energy) has units 1/nm^2 and
integrates to 1 over the whole image plane, exactly like the unit-sum PSF of
`PSFSimulator` (which is I * pixel_area). Metrics then follow the definitions
in `psfsim.metrics`:

    strehl             max_r I(r, z) / I_nomask(0, 0)   (peak, as PSFResult.strehl)
    strehl_axial       I(0, z) / I_nomask(0, 0)          (on axis)
    power_ropt         P(r_opt) = Integral_0^r_opt I 2 pi r dr
    total_power        P(total_radius)
    encircled_energy   P(r_opt) / P(total_radius)
    ee50_nm            radius holding 50 % of the *total* energy (=1)
    second_moment_nm2  <r^2> inside r_opt
    power_weighted_m2  second_moment * P_ref / P(r_opt), P_ref = no mask, z = 0

Depth budgets (mask independent)
--------------------------------
Parseval in z (at fixed r) gives, for ANY radially symmetric phase,

    Integral S_axial dz = (lambda/n) 2 (n/NA)^2 Int A^2 c du / (Int A du)^2
    Integral I(r,z) dz  = 2 pi (n/lambda) Int A^2 c J0^2 du / Int A^2 du

so the z-integral of the on-axis Strehl and of P(r_opt) cannot be changed by
the mask, only redistributed along z. `RadialPSF.budgets()` evaluates both.
Note: the *peak* Strehl and the *ratio* EE = P(r_opt)/P(10 um) are not
strictly invariant (the peak can leave the axis; the denominator loses light
beyond 10 um at large |z|), so they are reported but not used as the budget.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional, Sequence

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import j0

from .config import OpticalSystem
from .masks import build_mask
from .pupil import cos_theta, defocus_phase, pupil_amplitude, radial_pupil_grid
from .zernike import ModeSpec, resolve_mode, zernike_map

PhaseSpec = Optional[np.ndarray | Callable[[np.ndarray], np.ndarray]]

METRIC_KEYS = (
    "strehl",
    "strehl_axial",
    "ee50_nm",
    "encircled_energy",
    "second_moment_nm2",
    "power_weighted_m2",
    "power_ropt",
    "total_power",
)


@dataclass
class RadialStack:
    """Per-plane metrics of one mask over a depth sweep (radial model)."""

    z_nm: np.ndarray
    r_nm: np.ndarray
    intensity: np.ndarray  # (n_z, n_r), 1/nm^2, unit total power
    metrics: dict[str, np.ndarray]
    r_opt_nm: float
    total_radius_nm: float
    label: str = ""

    def __getitem__(self, key: str) -> np.ndarray:
        return self.metrics[key]

    def as_dict(self) -> dict[str, np.ndarray]:
        return dict(self.metrics)


@dataclass
class Budgets:
    """Mask-independent z-integrals (infinite range), in nm."""

    strehl_axial_nm: float  # Integral S_axial(z) dz
    power_ropt_nm: float  # Integral P(r_opt, z) dz
    r_opt_nm: float
    extra: dict = field(default_factory=dict)


class RadialPSF:
    """
    1D Debye evaluator for radially symmetric masks.

    Parameters
    ----------
    system : optical system (shared with PSFSimulator).
    n_pupil : number of equal-area pupil rings (midpoint rule in u = rho^2).
        2048 is converged to <1e-4 for |z| <= 8 um and 10 um radius.
    r_opt_nm, total_radius_nm : the EE radii (see psfsim.metrics).
    dr_nm : radial image sampling; adjusted so r_opt and total_radius are nodes.
    exact_defocus : exact (cos theta) or paraxial depth phase.
    """

    def __init__(
        self,
        system: OpticalSystem | None = None,
        n_pupil: int = 2048,
        r_opt_nm: float = 400.0,
        total_radius_nm: float = 10_000.0,
        dr_nm: float = 10.0,
        exact_defocus: bool = True,
    ) -> None:
        self.system = system or OpticalSystem()
        s = self.system
        self.n_pupil = int(n_pupil)
        self.r_opt_nm = float(r_opt_nm)
        self.total_radius_nm = float(total_radius_nm)

        # Pupil nodes: midpoints in u = rho^2 (equal area per node).
        self.du = 1.0 / self.n_pupil
        self.u = (np.arange(self.n_pupil) + 0.5) * self.du
        self.rho = np.sqrt(self.u)
        self.grid = radial_pupil_grid(self.rho)
        self.c = cos_theta(self.grid, s)
        self.amplitude = pupil_amplitude(self.grid, s)
        self.weights = self.amplitude * self.du / 2.0  # A rho drho
        self.energy = float(np.sum(self.amplitude**2) * self.du / 2.0)
        # Defocus is linear in z: evaluate once per nm.
        self.defocus_per_nm = defocus_phase(self.grid, s, 1000.0, exact=exact_defocus) / 1000.0

        # Image radii: r_opt and total_radius are exact nodes.
        n_opt = max(2, int(round(self.r_opt_nm / dr_nm)))
        n_opt += n_opt % 2  # even number of intervals for Simpson
        self.dr_nm = self.r_opt_nm / n_opt
        self._i_opt = n_opt
        self._i_tot = int(round(self.total_radius_nm / self.dr_nm))
        self.r_nm = np.arange(self._i_tot + 1) * self.dr_nm
        self._kr = 2.0 * np.pi * s.na / s.wavelength_nm
        self._J = j0(self._kr * self.rho[:, None] * self.r_nm[None, :])  # (n_pupil, n_r)

        self._norm = 2.0 * np.pi * (s.na / s.wavelength_nm) ** 2 / self.energy
        self.reference_peak = self._norm * float(np.sum(self.weights)) ** 2  # I_nomask(0,0)
        ref, _ = self._metrics(np.array([0.0]), np.zeros(self.n_pupil), ref_power=None)
        self.reference_power = float(ref["power_ropt"][0])

    # ------------------------------------------------------------------ #
    # Phase assembly
    # ------------------------------------------------------------------ #
    def mask_phase(self, mask: str | None, mask_params: Mapping | None = None) -> np.ndarray:
        """Evaluate a registered mask on the radial nodes (waves). Raises if
        the mask is not radially symmetric."""
        if mask is None or mask == "none":
            return np.zeros(self.n_pupil)
        phase = build_mask(mask, self.grid, self.system, mask_params)
        rot = radial_pupil_grid(self.rho)
        t = 1.2345
        rot.x, rot.y, rot.theta = self.rho * np.cos(t), self.rho * np.sin(t), np.full_like(self.rho, t)
        phase_rot = build_mask(mask, rot, self.system, mask_params)
        if not np.allclose(phase, phase_rot, atol=1e-9):
            raise ValueError(f"mask {mask!r} is not radially symmetric; use PSFSimulator")
        return np.asarray(phase, dtype=float)

    def aberration_phase(self, aberrations: Mapping[ModeSpec, float] | None) -> np.ndarray:
        """Radially symmetric (m = 0) Zernike aberrations, Noll-normalised RMS waves."""
        out = np.zeros(self.n_pupil)
        for mode, coeff in (aberrations or {}).items():
            n, m = resolve_mode(mode)
            if m != 0:
                raise ValueError(f"mode {mode} (n={n}, m={m}) is not radially symmetric")
            out += coeff * zernike_map((n, 0), self.rho, self.grid.theta, None, normalize=True)
        return out

    def total_phase(
        self,
        mask: str | None = None,
        mask_params: Mapping | None = None,
        aberrations: Mapping[ModeSpec, float] | None = None,
        phase_waves: PhaseSpec = None,
    ) -> np.ndarray:
        """Depth-independent pupil phase (waves) on the radial nodes. `phase_waves`
        may be an array on `self.rho` or a callable of rho, for free-form masks."""
        out = self.mask_phase(mask, mask_params) + self.aberration_phase(aberrations)
        if phase_waves is not None:
            out = out + (phase_waves(self.rho) if callable(phase_waves) else np.asarray(phase_waves))
        return out

    # ------------------------------------------------------------------ #
    # Core evaluation
    # ------------------------------------------------------------------ #
    def _pupil_fields(self, z_nm: np.ndarray, phase: np.ndarray) -> np.ndarray:
        """(n_z, n_pupil) complex weights w * exp(2 pi i psi(z))."""
        psi = phase[None, :] + z_nm[:, None] * self.defocus_per_nm[None, :]
        return self.weights[None, :] * np.exp(2j * np.pi * psi)

    def intensity(self, z_nm: Sequence[float], phase: np.ndarray, chunk: int = 256) -> np.ndarray:
        """I(r, z) on (z_nm, self.r_nm), 1/nm^2, unit total power per plane."""
        z = np.atleast_1d(np.asarray(z_nm, dtype=float))
        out = np.empty((z.size, self.r_nm.size))
        for i in range(0, z.size, chunk):
            e = self._pupil_fields(z[i : i + chunk], phase)
            f_re = e.real @ self._J
            f_im = e.imag @ self._J
            out[i : i + chunk] = self._norm * (f_re**2 + f_im**2)
        return out

    def axial_strehl(self, z_nm: Sequence[float], phase: np.ndarray) -> np.ndarray:
        """On-axis Strehl only (cheap; for optimizer inner loops)."""
        z = np.atleast_1d(np.asarray(z_nm, dtype=float))
        f = self._pupil_fields(z, phase).sum(axis=1)
        return np.abs(f) ** 2 / float(np.sum(self.weights)) ** 2

    def _metrics(self, z: np.ndarray, phase: np.ndarray, ref_power: float | None, intensity=None):
        inten = self.intensity(z, phase) if intensity is None else intensity
        r = self.r_nm
        cum = cumulative_simpson(inten * (2.0 * np.pi * r)[None, :], x=r, axis=1, initial=0.0)
        p_opt = cum[:, self._i_opt]
        p_tot = cum[:, self._i_tot]
        sl = slice(0, self._i_opt + 1)
        m2 = simpson(inten[:, sl] * (2.0 * np.pi * r[sl] ** 3)[None, :], x=r[sl], axis=1) / p_opt
        ee50 = np.full(z.size, np.nan)
        for k in range(z.size):
            j = int(np.searchsorted(cum[k], 0.5))
            if 0 < j < r.size:
                c0, c1 = cum[k, j - 1], cum[k, j]
                ee50[k] = r[j - 1] + (0.5 - c0) / (c1 - c0) * self.dr_nm
        return {
            "strehl": inten.max(axis=1) / self.reference_peak,
            "strehl_axial": inten[:, 0] / self.reference_peak,
            "ee50_nm": ee50,
            "encircled_energy": p_opt / p_tot,
            "second_moment_nm2": m2,
            "power_weighted_m2": m2 * ref_power / p_opt if ref_power else np.full(z.size, np.nan),
            "power_ropt": p_opt,
            "total_power": p_tot,
        }, inten

    def simulate(
        self,
        z_nm: Sequence[float],
        mask: str | None = None,
        mask_params: Mapping | None = None,
        aberrations: Mapping[ModeSpec, float] | None = None,
        phase_waves: PhaseSpec = None,
        label: str = "",
    ) -> RadialStack:
        """All metrics over a depth sweep. Same mask arguments as PSFSimulator."""
        z = np.atleast_1d(np.asarray(z_nm, dtype=float))
        phase = self.total_phase(mask, mask_params, aberrations, phase_waves)
        metrics, inten = self._metrics(z, phase, self.reference_power)
        return RadialStack(
            z_nm=z,
            r_nm=self.r_nm,
            intensity=inten,
            metrics=metrics,
            r_opt_nm=self.r_opt_nm,
            total_radius_nm=self.total_radius_nm,
            label=label or (mask or "none"),
        )

    # ------------------------------------------------------------------ #
    # Budgets
    # ------------------------------------------------------------------ #
    def budgets(self) -> Budgets:
        """Infinite-range z-integrals of S_axial and P(r_opt). Mask independent
        by construction (no phase enters); see the module docstring."""
        s = self.system
        a, c, du = self.amplitude, self.c, self.du
        b_s = (s.wavelength_nm / s.n_medium) * 2.0 * (s.n_medium / s.na) ** 2 * (
            np.sum(a**2 * c) * du / (np.sum(a) * du) ** 2
        )
        sl = slice(0, self._i_opt + 1)
        r = self.r_nm[sl]
        jj = self._J[:, sl] ** 2
        iz = 2.0 * np.pi * (s.n_medium / s.wavelength_nm) * ((a**2 * c) @ jj) / np.sum(a**2)
        b_p = simpson(iz * 2.0 * np.pi * r, x=r)
        paraxial = 2.0 * s.n_medium * s.wavelength_nm / s.na**2
        uniform_g = s.wavelength_nm / (s.n_medium * (1.0 - s.cos_min))
        return Budgets(
            strehl_axial_nm=float(b_s),
            power_ropt_nm=float(b_p),
            r_opt_nm=self.r_opt_nm,
            extra={"paraxial_strehl_nm": paraxial, "uniform_g_strehl_nm": uniform_g},
        )
