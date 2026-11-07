"""
Phase masks displayed on the SLM.

Every mask is a function `(grid, system, **params) -> phase map in waves`,
registered by name so that PSF-engineering work later only has to add a new
function, not touch the pipeline:

    @register_mask("my_mask")
    def _my_mask(grid, system, strength_waves=1.0):
        return strength_waves * something(grid.rho, grid.theta)

Masks return *ideal, continuous, unwrapped* phase. Wrapping, pixelation, and
quantisation are the SLM's business and live in `psfsim.slm`.

The base case requested — primary spherical — is `"spherical"`.
"""

from __future__ import annotations

from typing import Callable, Dict, Mapping

import numpy as np

from .config import OpticalSystem
from .pupil import PupilGrid
from .zernike import ModeSpec, zernike_map

MaskFunction = Callable[..., np.ndarray]
MASK_REGISTRY: Dict[str, MaskFunction] = {}


def register_mask(name: str) -> Callable[[MaskFunction], MaskFunction]:
    """Decorator registering a mask under `name`."""

    def wrapper(fn: MaskFunction) -> MaskFunction:
        if name in MASK_REGISTRY:
            raise KeyError(f"mask {name!r} is already registered")
        MASK_REGISTRY[name] = fn
        return fn

    return wrapper


def build_mask(
    name: str | None,
    grid: PupilGrid,
    system: OpticalSystem,
    params: Mapping | None = None,
) -> np.ndarray:
    """Look up and evaluate a mask. `name=None` (or "none") gives a flat pupil."""
    if name is None or name == "none":
        return np.zeros_like(grid.rho)
    if name not in MASK_REGISTRY:
        raise KeyError(
            f"unknown mask {name!r}; available: {sorted(MASK_REGISTRY)}"
        )
    phase = MASK_REGISTRY[name](grid, system, **(params or {}))
    return np.where(grid.mask, phase, 0.0)


def available_masks() -> list[str]:
    return sorted(MASK_REGISTRY)


# --------------------------------------------------------------------------- #
# Masks
# --------------------------------------------------------------------------- #
@register_mask("spherical")
def _spherical(
    grid: PupilGrid,
    system: OpticalSystem,
    strength_waves: float = 1.5,
    normalized: bool = True,
) -> np.ndarray:
    """
    Primary spherical aberration, Z(n=4, m=0) — the base EDOF mask.

    Deliberately adding spherical aberration trades peak intensity for axial
    invariance: the focus is smeared into a tube instead of a point, so the PSF
    stays usable over several times the native DOF. It is the cheapest EDOF mask
    to display (radially symmetric, no vortex) and its strength is a single
    scalar, which makes it a good baseline before moving to axial-encoding PSFs.

    strength_waves : RMS wavefront in waves if `normalized` (Noll-normalised
        Z11, so RMS = strength_waves), otherwise the plain 6r^4-6r^2+1 amplitude.
        Typical useful range is ~0.5-3 waves RMS; past that the Strehl drop
        starts costing more SNR than the depth gain buys.
    """
    return strength_waves * zernike_map(
        (4, 0), grid.rho, grid.theta, grid.mask, normalize=normalized
    )


@register_mask("annular")
def _annular(
    grid: PupilGrid,
    system: OpticalSystem,
    radii: tuple[float, ...] = (0.76, 0.93),
    phases_waves: tuple[float, ...] = (0.0, 0.5, 0.0),
) -> np.ndarray:
    """
    Concentric zones of constant phase (binary / multi-level annular mask).

    `radii` are the zone boundaries in normalised pupil radius (strictly
    increasing, inside (0, 1)); `phases_waves` gives one phase per zone, so it
    has len(radii) + 1 entries. Zone k covers radii[k-1] <= rho < radii[k].

    Defaults are the three-zone binary 0 / pi / 0 mask with boundaries at
    rho = 0.76 and 0.93. The pi-shifted ring interferes destructively with the
    central disk at nominal focus and less so away from it, redistributing
    energy along z. Practical advantages over a continuous mask: only two gray
    levels, and no phase gradient for the SLM to resolve except at the two
    edges. 0.5 waves is exactly representable on an 8-bit SLM (level 128), so
    there is no quantisation error at the default bit depth.
    """
    radii = tuple(float(r) for r in radii)
    phases = tuple(float(p) for p in phases_waves)
    if len(phases) != len(radii) + 1:
        raise ValueError(f"need len(radii) + 1 = {len(radii) + 1} phases, got {len(phases)}")
    if any(b <= a for a, b in zip(radii, radii[1:])) or not all(0 < r < 1 for r in radii):
        raise ValueError("radii must be strictly increasing and inside (0, 1)")
    zone = np.searchsorted(np.asarray(radii), grid.rho, side="right")
    return np.asarray(phases)[zone]


def equal_area_radii(n_zones: int) -> tuple[float, ...]:
    """Boundaries splitting the unit pupil into `n_zones` rings of equal area,
    r_k = sqrt(k / n). Equal area means equal transmitted energy per zone for a
    uniform pupil. Two zones gives a single boundary at 1/sqrt(2) ~ 0.7071."""
    if n_zones < 1:
        raise ValueError("need at least one zone")
    return tuple(float(np.sqrt(k / n_zones)) for k in range(1, n_zones))


@register_mask("zoned_zernike")
def _zoned_zernike(
    grid: PupilGrid,
    system: OpticalSystem,
    radii: tuple[float, ...] = (1 / np.sqrt(2),),
    zone_coefficients: tuple[Mapping[ModeSpec, float], ...] = ({4: 0.5}, {4: -0.5}),
    normalize: bool = False,
    local_coords: bool = False,
) -> np.ndarray:
    """
    Layer-cake mask whose zones each carry their own Zernike phase profile,
    instead of the constant phase of the `annular` mask.

    radii : zone boundaries in normalised pupil radius, strictly increasing.
    zone_coefficients : one {mode: coefficient in waves} dict per zone
        (len(radii) + 1 of them), innermost first.
    normalize : False (default) uses the plain polynomials, so a coefficient is
        the *peak amplitude* of the term — e.g. defocus is c*(2 rho^2 - 1),
        spanning -c..+c over the full pupil. True uses Noll normalisation
        (coefficient = RMS over the full pupil), as the other masks do.
    local_coords : False (default) evaluates every zone's polynomial in the
        global pupil coordinate rho, so each zone is literally a piece of the
        full-pupil Zernike. True re-maps each zone onto its own unit disk by
        an area-preserving change of radius, rho' = sqrt((rho^2 - a^2) /
        (b^2 - a^2)), so every zone carries a complete copy of the mode.

    Default: two equal-area zones, inner +0.5 waves (pi) and outer -0.5 waves
    (-pi) of defocus. Because the boundary sits at rho^2 = 1/2, exactly where
    2 rho^2 - 1 crosses zero, this phase is continuous: it is -pi * |2 rho^2 - 1|,
    i.e. the two halves of the pupil receive equal and opposite defocus. Each
    half then focuses at a different depth, giving a bifocal PSF with half the
    energy in each focus.
    """
    radii = tuple(float(r) for r in radii)
    if len(zone_coefficients) != len(radii) + 1:
        raise ValueError(
            f"need len(radii) + 1 = {len(radii) + 1} coefficient dicts, "
            f"got {len(zone_coefficients)}"
        )
    if any(b <= a for a, b in zip(radii, radii[1:])) or not all(0 < r < 1 for r in radii):
        raise ValueError("radii must be strictly increasing and inside (0, 1)")

    edges = (0.0,) + radii + (1.0,)
    zone = np.searchsorted(np.asarray(radii), grid.rho, side="right")
    out = np.zeros_like(grid.rho)
    for k, coeffs in enumerate(zone_coefficients):
        sel = (zone == k) & grid.mask
        if not coeffs or not sel.any():
            continue
        if local_coords:
            a, b = edges[k], edges[k + 1]
            rho = np.sqrt(np.clip((grid.rho**2 - a**2) / (b**2 - a**2), 0.0, 1.0))
        else:
            rho = grid.rho
        for mode, c in coeffs.items():
            z = zernike_map(mode, rho, grid.theta, None, normalize=normalize)
            out[sel] += c * z[sel]
    return out


@register_mask("zernike")
def _zernike_combo(
    grid: PupilGrid,
    system: OpticalSystem,
    coefficients: Mapping[ModeSpec, float] | None = None,
) -> np.ndarray:
    """
    Arbitrary Zernike combination, e.g. {11: 1.5, 22: -0.4}.

    This is the general-purpose entry point for both AO correction patterns and
    tuned EDOF masks (the "Zernike tuning" family, where a few radial orders are
    optimised jointly for axial invariance).
    """
    out = np.zeros_like(grid.rho)
    for mode, coeff in (coefficients or {}).items():
        out += coeff * zernike_map(mode, grid.rho, grid.theta, grid.mask)
    return out


@register_mask("defocus_sweep")
def _defocus_sweep(
    grid: PupilGrid,
    system: OpticalSystem,
    depth_range_nm: float = 2000.0,
) -> np.ndarray:
    """
    Cubic-in-rho^2 phase equivalent to sweeping focus across `depth_range_nm`
    during one exposure (the "SPED"-style focus-sweep EDOF). Included as a
    second baseline because it is the mask whose depth response is easiest to
    predict analytically: the axial MTF is roughly the integral of the in-focus
    OTF over the sweep.
    """
    w_max = (system.na**2 * depth_range_nm) / (2 * system.n_medium * system.wavelength_nm)
    return np.where(grid.mask, w_max * grid.rho**4 / 3.0, 0.0)


@register_mask("cubic")
def _cubic(
    grid: PupilGrid,
    system: OpticalSystem,
    alpha_waves: float = 10.0,
) -> np.ndarray:
    """
    Classic cubic (wavefront-coding) EDOF mask, alpha*(x^3 + y^3) in normalised
    pupil coordinates. Non-symmetric PSF that is highly defocus-invariant but
    requires deconvolution — kept here for comparison against the spherical
    baseline.
    """
    return alpha_waves * (grid.x**3 + grid.y**3)


@register_mask("tilt")
def _tilt(
    grid: PupilGrid,
    system: OpticalSystem,
    x_waves: float = 0.0,
    y_waves: float = 0.0,
) -> np.ndarray:
    """
    Linear ramp. Not an imaging mask: it is the carrier you add to push the
    modulated light off-axis and dump the SLM's unmodulated zero order, so it
    belongs in the same pipeline stage as the mask itself.
    """
    return x_waves * grid.x + y_waves * grid.y
