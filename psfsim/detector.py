"""
Camera sampling of the fine-grid PSF.

NOTE (Phase 0): reconstructed; the original module was not in the skeleton.
Camera pixels integrate the fine PSF over their area (exact overlap of square
pixels, no noise). Noise models come with Phase 5.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CameraImage:
    image: np.ndarray
    pixel_size_nm: float

    @property
    def shape(self) -> tuple[int, ...]:
        return self.image.shape

    @property
    def x_nm(self) -> np.ndarray:
        n = self.image.shape[-1]
        return (np.arange(n) - n // 2) * self.pixel_size_nm


def _overlap_matrix(n_fine: int, fine_nm: float, n_cam: int, cam_nm: float) -> np.ndarray:
    """W[i, j] = fraction of fine pixel j that lies inside camera pixel i."""
    f_lo = (np.arange(n_fine) - n_fine // 2 - 0.5) * fine_nm
    f_hi = f_lo + fine_nm
    c_lo = (np.arange(n_cam) - n_cam // 2 - 0.5) * cam_nm
    c_hi = c_lo + cam_nm
    ov = np.minimum(c_hi[:, None], f_hi[None, :]) - np.maximum(c_lo[:, None], f_lo[None, :])
    return np.clip(ov, 0.0, None) / fine_nm


def downsample_to_camera(
    psf: np.ndarray, sim_pixel_nm: float, camera_pixel_nm: float, fov_nm: float | None = None
) -> CameraImage:
    n_fine = psf.shape[0]
    extent = n_fine * sim_pixel_nm if fov_nm is None else min(fov_nm, n_fine * sim_pixel_nm)
    n_cam = int(np.floor(extent / camera_pixel_nm))
    n_cam -= (n_cam + 1) % 2  # odd, so the optical axis is a pixel centre
    n_cam = max(n_cam, 1)
    w = _overlap_matrix(n_fine, sim_pixel_nm, n_cam, camera_pixel_nm)
    keep = w.any(axis=0)
    wk = w[:, keep]
    img = wk @ psf[np.ix_(keep, keep)] @ wk.T
    return CameraImage(image=img, pixel_size_nm=camera_pixel_nm)


def crop_center(img: np.ndarray, pixel_nm: float, fov_nm: float) -> tuple[np.ndarray, np.ndarray]:
    n = img.shape[0]
    half = min(int(round(fov_nm / pixel_nm / 2)), n // 2 - 1 + n % 2)
    c = n // 2
    sl = slice(c - half, c + half + 1)
    return img[sl, sl], (np.arange(-half, half + 1)) * pixel_nm
