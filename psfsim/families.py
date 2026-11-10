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

import re
from dataclasses import dataclass, field

import numpy as np

from .config import OpticalSystem
from .masks import equal_dc_radii, t_to_rho

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


def family_by_name(name: str, system: OpticalSystem | None = None) -> Family:
    """Rebuild a family from the name it is logged under (e.g. in
    docs/phase2/results.json): "zoned_focus_4_sweep", "annular_binary_3", ..."""
    fixed = {"spherical+defocus": SphericalDefocus, "zernike_m0": ZernikeM0,
             "zoned_zernike2": ZonedZernike2}
    if name in fixed:
        return fixed[name](system)
    m = re.fullmatch(r"(annular_binary|layer_cake|zoned_focus)_(\d+)(_sweep)?", name)
    if not m or (m.group(3) and m.group(1) != "zoned_focus"):
        raise KeyError(f"unknown family {name!r}")
    n = int(m.group(2))
    if m.group(1) == "zoned_focus":
        return ZonedFocus(n, sweep=bool(m.group(3)), system=system)
    return {"annular_binary": AnnularBinary, "layer_cake": LayerCake}[m.group(1)](n, system)
