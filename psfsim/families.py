"""
Parameterised mask families for optimisation (Phase 1-2).

Each family maps a parameter vector x (box-bounded by `lo`, `hi`) to a
registered mask and its parameters, so any optimiser (grid, CMA-ES, NSGA-II,
Bayesian) can drive it through `Evaluator.evaluate_mask(*family.mask(x))`.

Zone edges are parameterised by positive widths in t = (1 - cos theta) /
(1 - cos theta_max) (equal t = equal depth of focus per zone), normalised to
sum to one. This keeps the edges ordered for any x without constraints; the
overall scale of the widths is a redundant direction, which CMA-ES tolerates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import OpticalSystem
from .masks import equal_dc_radii, rho_to_t, t_to_rho
from .zernike import radial_polynomial

DEPTH_BOUND_NM = 4000.0
SWEEP_BOUND_NM = 4000.0
MIN_WIDTH = 0.05  # relative zone width lower bound (before normalisation)


def widths_to_radii(w, system: OpticalSystem) -> tuple[float, ...]:
    w = np.asarray(w, dtype=float)
    t = np.cumsum(w)[:-1] / w.sum()
    return tuple(float(r) for r in t_to_rho(t, system))


@dataclass
class Family:
    name: str
    names: list[str]
    lo: np.ndarray
    hi: np.ndarray
    system: OpticalSystem = field(default_factory=OpticalSystem)

    @property
    def dim(self) -> int:
        return len(self.names)

    def mask(self, x) -> tuple[str, dict]:
        raise NotImplementedError

    def describe(self, x) -> str:
        m, p = self.mask(x)

        def fmt(v):
            if isinstance(v, dict):
                return "{" + ", ".join(f"{k}: {fmt(u)}" for k, u in v.items()) + "}"
            if isinstance(v, (tuple, list, np.ndarray)):
                return "(" + ", ".join(fmt(u) for u in v) + ")"
            return f"{float(v):.4g}"

        return f"{m}: " + ", ".join(f"{k}={fmt(v)}" for k, v in p.items())


class SphericalDefocus(Family):
    """Z11 spherical + Z4 defocus (Noll RMS waves): the re-centred spherical baseline."""

    def __init__(self, system=None):
        super().__init__("spherical+defocus", ["z4", "z11"], np.array([-1.2, -1.0]),
                         np.array([1.2, 1.0]), system or OpticalSystem())

    def mask(self, x):
        return "zernike", {"coefficients": {4: float(x[0]), 11: float(x[1])}}


class ZernikeM0(Family):
    """Z4, Z11, Z22, Z37 (Noll RMS waves): Phase 1 family A."""

    def __init__(self, system=None):
        super().__init__("zernike_m0", ["z4", "z11", "z22", "z37"],
                         np.array([-1.2, -1.0, -0.6, -0.4]), np.array([1.2, 1.0, 0.6, 0.4]),
                         system or OpticalSystem())

    def mask(self, x):
        return "zernike", {"coefficients": dict(zip((4, 11, 22, 37), map(float, x)))}


class ZonedZernike2(Family):
    """Current zoned_zernike, two zones, free edge, paraxial defocus per zone,
    piston on the outer zone: Phase 1 family B."""

    def __init__(self, system=None):
        super().__init__("zoned_zernike2", ["rho1", "c_in", "c_out", "piston"],
                         np.array([0.3, -1.5, -1.5, -1.0]), np.array([0.95, 1.5, 1.5, 2.0]),
                         system or OpticalSystem())

    def mask(self, x):
        rho1, c_in, c_out, p = map(float, x)
        return "zoned_zernike", {"radii": (rho1,), "zone_coefficients": ({4: c_in}, {4: c_out, 1: p})}


class AnnularBinary(Family):
    """Binary 0 / pi annular mask with `n_edges` free edges (alternating phase)."""

    def __init__(self, n_edges: int, system=None):
        n = n_edges + 1
        super().__init__(f"annular_binary_{n_edges}", [f"w{k}" for k in range(n)],
                         np.full(n, MIN_WIDTH), np.ones(n), system or OpticalSystem())
        self.n_edges = n_edges

    def mask(self, x):
        radii = widths_to_radii(x, self.system)
        phases = tuple(0.5 * (k % 2) for k in range(self.n_edges + 1))
        return "annular", {"radii": radii, "phases_waves": phases}


class LayerCake(Family):
    """Coherent layer cake: N equal-dc zones, no defocus, free pistons (zone 0 = 0)."""

    def __init__(self, n_zones: int, system=None):
        super().__init__(f"layer_cake_{n_zones}", [f"p{k}" for k in range(1, n_zones)],
                         np.full(n_zones - 1, -1.0), np.full(n_zones - 1, 2.0), system or OpticalSystem())
        self.n_zones = n_zones
        self.radii = equal_dc_radii(n_zones, self.system)

    def mask(self, x):
        return "zoned_focus", {"radii": self.radii, "depths_nm": (0.0,) * self.n_zones,
                               "pistons_waves": (0.0,) + tuple(map(float, x))}


class ZonedFocus(Family):
    """Phase 2 family: N zones with free edges (widths in t), exact-defocus target
    depths z_k, pistons p_k (zone 0 fixed at 0) and optional sweeps s_k."""

    def __init__(self, n_zones: int, sweep: bool = False, system=None):
        n = n_zones
        names = [f"w{k}" for k in range(n)] + [f"z{k}" for k in range(n)] + [f"p{k}" for k in range(1, n)]
        lo = [MIN_WIDTH] * n + [-DEPTH_BOUND_NM] * n + [-1.0] * (n - 1)
        hi = [1.0] * n + [DEPTH_BOUND_NM] * n + [2.0] * (n - 1)
        if sweep:
            names += [f"s{k}" for k in range(n)]
            lo += [-SWEEP_BOUND_NM] * n
            hi += [SWEEP_BOUND_NM] * n
        super().__init__(f"zoned_focus_{n}{'_sweep' if sweep else ''}", names, np.array(lo, float),
                         np.array(hi, float), system or OpticalSystem())
        self.n_zones, self.sweep = n, sweep

    def mask(self, x):
        n = self.n_zones
        x = np.asarray(x, dtype=float)
        p = {"radii": widths_to_radii(x[:n], self.system),
             "depths_nm": tuple(x[n:2 * n]),
             "pistons_waves": (0.0,) + tuple(x[2 * n:3 * n - 1])}
        if self.sweep:
            p["sweeps_nm"] = tuple(x[3 * n - 1:4 * n - 1])
        return "zoned_focus", p


def noll_m0(n: int) -> int:
    """Noll index of the radially symmetric Zernike Z(n, 0): 1, 4, 11, 22, 37, ..."""
    return n * (n + 1) // 2 + 1


class ZonedZernikeRadial(Family):
    """N zones with free edges (widths in t); every zone carries its own set of
    radially symmetric Zernikes Z(n, 0), n = 0, 2, ..., `max_order`, i.e. piston,
    defocus, primary spherical and (max_order = 6) secondary spherical.

    The polynomials live in the zone's own unit-disk coordinate
    rho' = sqrt((rho^2 - a^2) / (b^2 - a^2)) (`zoned_zernike`, local_coords) and are
    not normalised, so a coefficient is the peak amplitude of that term inside the
    zone, in waves. rho'^2 is affine in rho^2, so the span is the same as with
    full-pupil Zernikes cut into zones (any even polynomial of degree <= max_order
    per zone); the local form is used because the terms are orthogonal over the
    zone, which is what the optimiser needs. Zone 0's piston is the global piston
    and is fixed at 0.

    `max_order = 2` is the paraxial "zoned defocus + piston" mask. With
    `max_order = 6` the family contains `zoned_focus` (exact defocus, with or
    without sweep) to ~1e-3 waves: see `from_phase`.
    """

    COEF_BOUND = {0: (-1.0, 2.0), 2: (-2.0, 2.0), 4: (-1.0, 1.0), 6: (-0.5, 0.5), 8: (-0.25, 0.25)}

    def __init__(self, n_zones: int, max_order: int = 6, system=None):
        if max_order % 2 or not 2 <= max_order <= 8:
            raise ValueError("max_order must be 2, 4, 6 or 8")
        n = n_zones
        self.orders = tuple(range(0, max_order + 1, 2))
        self.terms = [(k, o) for k in range(n) for o in self.orders if (k, o) != (0, 0)]
        names = [f"w{k}" for k in range(n)] + [f"c{o}_{k}" for k, o in self.terms]
        lo = [MIN_WIDTH] * n + [self.COEF_BOUND[o][0] for _, o in self.terms]
        hi = [1.0] * n + [self.COEF_BOUND[o][1] for _, o in self.terms]
        super().__init__(f"zoned_zernike{max_order}_{n}", names, np.array(lo, float), np.array(hi, float),
                         system or OpticalSystem())
        self.n_zones, self.max_order = n, max_order

    def mask(self, x):
        n = self.n_zones
        x = np.asarray(x, dtype=float)
        coeffs = [dict() for _ in range(n)]
        for (k, o), v in zip(self.terms, x[n:]):
            coeffs[k][noll_m0(o)] = float(v)
        return "zoned_zernike", {"radii": widths_to_radii(x[:n], self.system),
                                 "zone_coefficients": tuple(coeffs), "local_coords": True}

    def from_phase(self, radii, rho: np.ndarray, phase_waves: np.ndarray) -> np.ndarray:
        """Parameter vector whose mask is the least-squares fit, zone by zone, of a
        radial phase profile `phase_waves` sampled on `rho` (equal-area nodes, e.g.
        `RadialPSF.rho`) with zone edges `radii`. Pistons are taken relative to
        zone 0 and wrapped into [0, 1) wave. Used to warm-start from another family."""
        n = self.n_zones
        radii = np.asarray(radii, dtype=float)
        if radii.size != n - 1:
            raise ValueError(f"need {n - 1} radii")
        edges = np.concatenate([[0.0], radii, [1.0]])
        fit = {}
        for k in range(n):
            a, b = edges[k], edges[k + 1]
            sel = (rho >= a) & (rho < b)
            r = np.sqrt(np.clip((rho[sel] ** 2 - a * a) / (b * b - a * a), 0.0, 1.0))
            basis = np.stack([radial_polynomial(o, 0, r) for o in self.orders], axis=1)
            co = np.linalg.lstsq(basis, np.asarray(phase_waves)[sel], rcond=None)[0]
            fit.update({(k, o): c for o, c in zip(self.orders, co)})
        p0 = fit[(0, 0)]
        t = np.concatenate([[0.0], rho_to_t(radii, self.system), [1.0]])
        w = np.diff(t)
        vals = [np.mod(fit[(k, o)] - p0, 1.0) if o == 0 else fit[(k, o)] for k, o in self.terms]
        return np.clip(np.concatenate([w / w.max(), vals]), self.lo, self.hi)
