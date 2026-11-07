"""
Tracking-oriented PSF metrics, evaluated per depth plane.

All metrics take a fine-grid PSF normalised to unit total energy (what
`PSFSimulator.simulate` returns in `PSFResult.psf`) and an "optimal PSF radius"
r_opt: the largest spot that a single-particle-tracking pipeline tolerates
before localisation precision degrades too much.

Centre
------
Both r_opt metrics are taken about the PSF's intensity-weighted centroid,
computed over the `total_radius_nm` disk around the optical axis. For the
radially symmetric masks studied so far this is the optical axis, but using the
centroid keeps the metrics honest for asymmetric (e.g. double-helix) PSFs
later. The brightest pixel is deliberately *not* used as the centre: far from
focus the peak of a defocused PSF often sits on a ring.

Sub-pixel accuracy
------------------
r_opt is typically only ~7 fine-grid samples, so a hard pixel-centre disk would
quantise the encircled energy in visible steps. Each pixel near the boundary is
therefore weighted by the fraction of its area inside the disk (supersampled),
and its r^2 is the mean over the covered part of the pixel.

Metrics
-------
encircled_energy     P(r_opt) / P(total_radius): fraction of the plane's light
                     that lands inside r_opt. Higher is better.
second_moment_nm2    <r^2> within r_opt, intensity weighted, about the centroid.
power_weighted_m2    second_moment * (P_ref / P(r_opt)), where P_ref is the
                     r_opt power of the no-mask PSF at focus. It penalises a
                     spot that is compact only because most light was thrown
                     outside r_opt. Lower is better; nm^2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass
class PlaneMetrics:
    centroid_nm: tuple[float, float]
    total_power: float  # energy within total_radius
    power_ropt: float  # energy within r_opt
    encircled_energy: float  # power_ropt / total_power
    second_moment_nm2: float  # <r^2> within r_opt
    power_weighted_m2: float  # second_moment * P_ref / power_ropt (nan if no ref)


def _coverage_window(
    shape: tuple[int, int],
    pixel_nm: float,
    center_px: tuple[float, float],
    radius_nm: float,
    supersample: int = 8,
):
    """
    For a disk of `radius_nm` about a (row, col) centre in pixel units, return
    (slice_r, slice_c, weight, mean_r2_nm2) over a small window: `weight` is
    the fraction of each pixel's area inside the disk, `mean_r2_nm2` the mean
    squared distance over the covered part of each pixel.
    """
    r_px = radius_nm / pixel_nm
    cy, cx = center_px
    r0 = max(int(np.floor(cy - r_px)) - 1, 0)
    r1 = min(int(np.ceil(cy + r_px)) + 2, shape[0])
    c0 = max(int(np.floor(cx - r_px)) - 1, 0)
    c1 = min(int(np.ceil(cx + r_px)) + 2, shape[1])

    s = supersample
    offs = (np.arange(s) + 0.5) / s - 0.5  # sub-sample offsets within a pixel
    rows = np.arange(r0, r1)[:, None] + offs[None, :]  # (nr, s)
    cols = np.arange(c0, c1)[:, None] + offs[None, :]  # (nc, s)
    dy2 = ((rows - cy) * pixel_nm) ** 2  # (nr, s)
    dx2 = ((cols - cx) * pixel_nm) ** 2  # (nc, s)
    # (nr, s, nc, s) squared distance of every sub-sample
    d2 = dy2[:, :, None, None] + dx2[None, None, :, :]
    inside = d2 <= radius_nm**2
    n_in = inside.sum(axis=(1, 3))
    weight = n_in / s**2
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_r2 = np.where(n_in > 0, (d2 * inside).sum(axis=(1, 3)) / np.maximum(n_in, 1), 0.0)
    return slice(r0, r1), slice(c0, c1), weight, mean_r2


def centroid_px(psf: np.ndarray, pixel_nm: float, radius_nm: float) -> tuple[float, float]:
    """Intensity-weighted centroid (row, col) within `radius_nm` of the grid centre."""
    n = psf.shape[0]
    c = n // 2
    idx = np.arange(n) - c
    yy, xx = np.meshgrid(idx, idx, indexing="ij")
    disk = (yy**2 + xx**2) * pixel_nm**2 <= radius_nm**2
    w = np.where(disk, psf, 0.0)
    tot = w.sum()
    if tot <= 0:
        return float(c), float(c)
    return float(c + (w * yy).sum() / tot), float(c + (w * xx).sum() / tot)


def plane_metrics(
    psf: np.ndarray,
    pixel_nm: float,
    r_opt_nm: float,
    total_radius_nm: float = 10_000.0,
    ref_power: float | None = None,
) -> PlaneMetrics:
    """All r_opt metrics for one PSF plane. See the module docstring."""
    n = psf.shape[0]
    c = n // 2
    cy, cx = centroid_px(psf, pixel_nm, total_radius_nm)

    # Total power: disk about the same centroid. Pixel-centre inclusion is
    # plenty at this scale (radius >> pixel).
    idx = np.arange(n)
    yy, xx = np.meshgrid(idx - cy, idx - cx, indexing="ij")
    total = float(psf[(yy**2 + xx**2) * pixel_nm**2 <= total_radius_nm**2].sum())

    sr, sc, w, r2 = _coverage_window(psf.shape, pixel_nm, (cy, cx), r_opt_nm)
    patch = psf[sr, sc]
    p_opt = float((patch * w).sum())
    m2 = float((patch * w * r2).sum() / p_opt) if p_opt > 0 else np.nan
    pwm2 = m2 * ref_power / p_opt if (ref_power is not None and p_opt > 0) else np.nan
    return PlaneMetrics(
        centroid_nm=((cy - c) * pixel_nm, (cx - c) * pixel_nm),
        total_power=total,
        power_ropt=p_opt,
        encircled_energy=p_opt / total if total > 0 else np.nan,
        second_moment_nm2=m2,
        power_weighted_m2=pwm2,
    )


def reference_power(sim, r_opt_nm: float, total_radius_nm: float = 10_000.0) -> float:
    """r_opt power of the no-mask, aberration-free PSF at the focal plane."""
    r = sim.simulate(z_nm=0.0)
    return plane_metrics(r.psf, sim.system.sim_pixel_nm, r_opt_nm, total_radius_nm).power_ropt


def depth_sweep_metrics(
    sim,
    z_values: Sequence[float],
    r_opt_nm: float,
    ref_power: float,
    total_radius_nm: float = 10_000.0,
    **simulate_kwargs,
) -> dict[str, np.ndarray]:
    """
    One mask over depth -> arrays of every metric. Single pass, no full results
    retained. `simulate_kwargs` go straight to `PSFSimulator.simulate`
    (mask, mask_params, aberrations, ...).
    """
    keys = ["strehl", "ee50_nm", "encircled_energy", "second_moment_nm2",
            "power_weighted_m2", "power_ropt", "total_power"]
    out = {k: [] for k in keys}
    px = sim.system.sim_pixel_nm
    for r in sim.iter_stack(z_values, **simulate_kwargs):
        m = plane_metrics(r.psf, px, r_opt_nm, total_radius_nm, ref_power)
        out["strehl"].append(r.strehl)
        out["ee50_nm"].append(r.encircled_energy_radius_nm)
        out["encircled_energy"].append(m.encircled_energy)
        out["second_moment_nm2"].append(m.second_moment_nm2)
        out["power_weighted_m2"].append(m.power_weighted_m2)
        out["power_ropt"].append(m.power_ropt)
        out["total_power"].append(m.total_power)
    return {k: np.asarray(v) for k, v in out.items()}


def contiguous_range(z: np.ndarray, ok: np.ndarray) -> tuple[float, float] | None:
    """
    Bounds of the contiguous run of `ok` that contains (or is nearest to) z = 0.
    Used for "usable depth range" summaries, where a disconnected region far from
    focus should not count.
    """
    ok = np.asarray(ok, dtype=bool)
    if not ok.any():
        return None
    i0 = int(np.argmin(np.abs(z)))
    if not ok[i0]:
        cand = np.flatnonzero(ok)
        i0 = int(cand[np.argmin(np.abs(z[cand]))])
    lo = hi = i0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(z) - 1 and ok[hi + 1]:
        hi += 1
    return float(z[lo]), float(z[hi])
