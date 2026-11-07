"""
SLM device model: pixelation, wrapping and quantisation of the mask phase.

NOTE (Phase 0): reconstructed and deliberately simple; the original module was
not in the skeleton. Phase 4 (hardware realism) is where this gets refined
(pixel crosstalk, LUT errors, fringing).

Steps, in order:
1. pixelate: the ideal phase is held constant over each SLM pixel, as imaged
   onto the pupil ("center" = value at the sample nearest the pixel centre,
   "mean" = average over the samples in the pixel). Skipped when the device
   is finer than the simulation grid.
2. wrap: modulo one wave (if system.slm_wrap).
3. quantise: 2**bit_depth levels per wave (if system.slm_bit_depth).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import OpticalSystem
from .pupil import PupilGrid


@dataclass
class SLMResult:
    applied_waves: np.ndarray
    samples_per_slm_pixel: float
    resolves_device_pixels: bool
    undersampled: bool
    max_step_waves_per_pixel: float


def _pixelate(phase: np.ndarray, grid: PupilGrid, system: OpticalSystem, mode: str) -> np.ndarray:
    pitch = 2.0 / system.slm_pixels_across_pupil  # SLM pixel in normalised pupil units
    ix = np.floor(grid.x / pitch).astype(np.int64)
    iy = np.floor(grid.y / pitch).astype(np.int64)
    ix -= ix.min()
    iy -= iy.min()
    label = (iy * (ix.max() + 1) + ix).ravel()
    inside = grid.mask.ravel()
    flat = phase.ravel()
    n_lab = int(label.max()) + 1
    if mode == "mean":
        s = np.bincount(label[inside], weights=flat[inside], minlength=n_lab)
        cnt = np.bincount(label[inside], minlength=n_lab)
        val = np.where(cnt > 0, s / np.maximum(cnt, 1), 0.0)
    elif mode == "center":
        cx = (np.floor(grid.x / pitch) + 0.5) * pitch
        cy = (np.floor(grid.y / pitch) + 0.5) * pitch
        d = np.hypot(grid.x - cx, grid.y - cy).ravel()
        order = np.lexsort((d, label))  # by label, then distance to centre
        first = np.ones(order.size, dtype=bool)
        first[1:] = label[order][1:] != label[order][:-1]
        val = np.zeros(n_lab)
        val[label[order][first]] = flat[order][first]
    else:
        raise ValueError("slm sampling mode must be 'center' or 'mean'")
    return np.where(grid.mask, val[label].reshape(phase.shape), 0.0)


def apply_slm(
    phase_waves: np.ndarray,
    grid: PupilGrid,
    system: OpticalSystem,
    mode: str = "center",
    pixelate_enabled: bool = True,
) -> SLMResult:
    samples = 2.0 * grid.pupil_radius_px / system.slm_pixels_across_pupil
    resolves = bool(np.isfinite(samples) and samples >= 1.0)

    # Largest phase change across one SLM pixel (sampling check).
    if phase_waves.ndim == 2:
        gy, gx = np.gradient(phase_waves)
        g = np.hypot(gx, gy)[grid.mask]
        step = float(g.max() * samples) if g.size and np.isfinite(samples) else 0.0
    else:
        step = 0.0

    out = phase_waves
    if pixelate_enabled and resolves and phase_waves.ndim == 2:
        out = _pixelate(out, grid, system, mode)
    if system.slm_wrap:
        out = np.mod(out, 1.0)
    if system.slm_bit_depth:
        levels = 2 ** int(system.slm_bit_depth)
        out = np.round(out * levels) / levels
        if system.slm_wrap:
            out = np.mod(out, 1.0)
    out = np.where(grid.mask, out, 0.0)
    return SLMResult(
        applied_waves=out,
        samples_per_slm_pixel=float(samples),
        resolves_device_pixels=resolves,
        undersampled=step > 0.5,
        max_step_waves_per_pixel=step,
    )
