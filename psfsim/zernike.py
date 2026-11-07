"""
Zernike polynomials on the unit pupil.

NOTE (Phase 0): reconstructed; the original module was not in the skeleton.

Modes are given either as a Noll index (int, 1-based) or as an (n, m) tuple.
m >= 0 is the cosine term, m < 0 the sine term. With `normalize=True`
(Noll normalisation) the RMS over the unit disk is 1, so a coefficient is an
RMS wavefront; with `normalize=False` the plain polynomial is used, e.g.
defocus is 2 rho^2 - 1 and primary spherical is 6 rho^4 - 6 rho^2 + 1.
"""

from __future__ import annotations

from math import factorial
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np

ModeSpec = Union[int, Tuple[int, int]]


def noll_to_nm(j: int) -> tuple[int, int]:
    """Noll index (1-based) -> (n, m)."""
    if j < 1:
        raise ValueError("Noll indices start at 1")
    n = 0
    j1 = j - 1
    while j1 > n:
        n += 1
        j1 -= n
    m = (-1) ** j * ((n % 2) + 2 * int((j1 + ((n + 1) % 2)) / 2))
    return n, int(m)


def resolve_mode(mode: ModeSpec) -> tuple[int, int]:
    if isinstance(mode, (int, np.integer)):
        return noll_to_nm(int(mode))
    n, m = mode
    n, m = int(n), int(m)
    if n < 0 or abs(m) > n or (n - abs(m)) % 2:
        raise ValueError(f"invalid Zernike mode (n={n}, m={m})")
    return n, m


def noll_range(j_max: int, j_min: int = 1) -> List[int]:
    return list(range(j_min, j_max + 1))


def radial_polynomial(n: int, m: int, rho: np.ndarray) -> np.ndarray:
    m = abs(m)
    out = np.zeros_like(rho, dtype=float)
    for k in range((n - m) // 2 + 1):
        coef = (-1) ** k * factorial(n - k) / (
            factorial(k) * factorial((n + m) // 2 - k) * factorial((n - m) // 2 - k)
        )
        out = out + coef * rho ** (n - 2 * k)
    return out


def zernike_map(
    mode: ModeSpec,
    rho: np.ndarray,
    theta: np.ndarray,
    mask: np.ndarray | None = None,
    normalize: bool = True,
) -> np.ndarray:
    n, m = resolve_mode(mode)
    z = radial_polynomial(n, m, rho)
    if m > 0:
        z = z * np.cos(m * theta)
    elif m < 0:
        z = z * np.sin(-m * theta)
    if normalize:
        z = z * (np.sqrt(n + 1) if m == 0 else np.sqrt(2 * (n + 1)))
    if mask is not None:
        z = np.where(mask, z, 0.0)
    return z


class ZernikeBasis:
    """Cached Zernike maps on a fixed grid, indexed by Noll index or (n, m)."""

    def __init__(
        self,
        modes: Iterable[ModeSpec],
        rho: np.ndarray,
        theta: np.ndarray,
        mask: np.ndarray | None,
    ) -> None:
        self._maps: Dict[tuple[int, int], np.ndarray] = {}
        for mode in modes:
            nm = resolve_mode(mode)
            self._maps[nm] = zernike_map(nm, rho, theta, mask, normalize=True)

    @property
    def modes(self) -> list[tuple[int, int]]:
        return list(self._maps)

    def __contains__(self, mode: ModeSpec) -> bool:
        return resolve_mode(mode) in self._maps

    def __getitem__(self, mode: ModeSpec) -> np.ndarray:
        return self._maps[resolve_mode(mode)]
