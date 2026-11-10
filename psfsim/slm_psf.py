"""
Pixel-level SLM model and PSF evaluation for lab-realistic mask comparison (Phase 4a).

The mask is drawn on the SLM pixel grid exactly as it is displayed:

* pixel grid: 2 R_mask x 2 R_mask pixels (R_mask = 362 for the 9.2 um SLM),
  optical axis on the corner between the four central pixels;
* the ideal (continuous) mask phase is evaluated at each pixel centre, in mask
  coordinates rho_mask = r / R_mask;
* the phase is "Fresnel-cut" into the SLM range [0, 2 pi) (mod one wave) and
  quantised to 2^bits levels (8-bit: level k -> 2 pi k / 256);
* optional pixel crosstalk ("fringing field" effect): the displayed phase map is
  blurred with a Gaussian of sigma `fringe_sigma_px` pixels. This is the one
  mechanism through which the 2 pi cut costs anything: an ideal cut is exact,
  because exp(i phi) is 2 pi periodic, but a smoothed 2 pi -> 0 edge spends ~1 px
  on wrong phases.

The true optical pupil (the back-focal-plane disk imaged onto the SLM) can have
a different radius `pupil_radius_px` than the drawn mask: 3.30 mm / 9.2 um =
358.7 px against a 362 px mask, i.e. the mask is ~0.9 % larger than the pupil.
The pupil, not the mask, sets NA, apodization and defocus.

Each SLM pixel is represented by `subsample` x `subsample` samples so that the
pixel staircase, the pupil edge and fringing are resolved. The field is
computed with a matrix Fourier transform (separable DFT kernels) directly on
two image grids:

* fine   : +-`fine_half_nm` at `fine_dx_nm`  -> Strehl, EE(r_opt), PW-M2, CRLB, images
* coarse : +-`coarse_half_nm` at `coarse_dx_nm` -> P(10 um) for the EE denominator, EE50

Intensities are normalised to unit total power by Parseval on the pupil, the
same convention as PSFSimulator / RadialPSF. With `ideal=True` the mask is
continuous (no pixels, no cut, no quantisation, no fringing), which reproduces
the radial model and serves as the reference.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np
from scipy.ndimage import gaussian_filter

from .config import OpticalSystem
from .crlb import Photometry
from .detector import _overlap_matrix
from .masks import build_mask
from .metrics import _coverage_window
from .objective import soft_min, trapezoid_weights
from .pupil import PupilGrid, defocus_phase, pupil_amplitude

SLM_PITCH_UM = 9.2
PUPIL_RADIUS_MM = 3.30


@dataclass(frozen=True)
class SLMDevice:
    mask_radius_px: float = 362.0  # drawn mask edge (rho = 1), SLM pixels
    pupil_radius_px: float = 362.0
    bit_depth: int | None = 8
    wrap: bool = True  # Fresnel cut into [0, 2 pi)
    fringe_sigma_px: float = 0.0
    ideal: bool = False  # continuous mask: ignore pixels, cut, quantisation, fringing
    label: str = ""

    def with_(self, **kw) -> "SLMDevice":
        return replace(self, **kw)


REAL_PUPIL_PX = PUPIL_RADIUS_MM * 1000.0 / SLM_PITCH_UM  # 358.7

DEVICES = {
    "ideal (continuous)": SLMDevice(ideal=True, label="ideal (continuous)"),
    "SLM, pupil = mask (362 px)": SLMDevice(label="SLM, pupil = mask (362 px)"),
    "SLM, pupil 3.30 mm (358.7 px)": SLMDevice(pupil_radius_px=REAL_PUPIL_PX,
                                               label="SLM, pupil 3.30 mm (358.7 px)"),
    "SLM, pupil 3.30 mm + fringing 0.5 px": SLMDevice(pupil_radius_px=REAL_PUPIL_PX, fringe_sigma_px=0.5,
                                                      label="SLM, pupil 3.30 mm + fringing 0.5 px"),
}


def _grid(xn: np.ndarray, yn: np.ndarray, scale: float) -> PupilGrid:
    """PupilGrid in units where `scale` pixels = 1."""
    x, y = xn / scale, yn / scale
    rho = np.hypot(x, y)
    return PupilGrid(x=x, y=y, rho=rho, theta=np.arctan2(y, x), mask=rho <= 1.0)


class SLMPSF:
    def __init__(self, system: OpticalSystem | None = None, device: SLMDevice | None = None,
                 subsample: int = 3, r_opt_nm: float = 400.0, total_radius_nm: float = 10_000.0,
                 fine_half_nm: float = 2500.0, fine_dx_nm: float = 25.0,
                 coarse_half_nm: float = 10_500.0, coarse_dx_nm: float = 100.0,
                 photometry: Photometry | None = None) -> None:
        self.system = system or OpticalSystem()
        self.device = device or SLMDevice()
        self.s = int(subsample)
        self.r_opt_nm, self.total_radius_nm = r_opt_nm, total_radius_nm
        d, s, sysm = self.device, self.s, self.system

        # sub-sample grid in SLM-pixel units, symmetric about the axis
        half_px = int(np.ceil(max(d.mask_radius_px, d.pupil_radius_px)))
        self.half_px = half_px
        u = (np.arange(2 * half_px * s) + 0.5) / s - half_px
        self.u_px = u
        X, Y = np.meshgrid(u, u, indexing="xy")
        self.pupil_grid = _grid(X, Y, d.pupil_radius_px)  # sets NA, apodization, defocus
        self.mask_grid = _grid(X, Y, d.mask_radius_px)  # where the mask is drawn
        # pixel-centre grid (for the displayed values)
        uc = np.arange(2 * half_px) + 0.5 - half_px
        Xc, Yc = np.meshgrid(uc, uc, indexing="xy")
        self.pixel_grid = _grid(Xc, Yc, d.mask_radius_px)

        self.aperture = self.pupil_grid.mask
        self.amplitude = pupil_amplitude(self.pupil_grid, sysm).astype(np.float32)
        self.defocus_per_nm = (defocus_phase(self.pupil_grid, sysm, 1000.0) / 1000.0).astype(np.float32)

        # pupil spatial frequencies (cycles / nm) and MFT kernels
        f = (sysm.na / sysm.wavelength_nm) * u / d.pupil_radius_px
        self.df = (sysm.na / sysm.wavelength_nm) / d.pupil_radius_px / s
        self.f = f
        self.x_fine = np.arange(-round(fine_half_nm / fine_dx_nm), round(fine_half_nm / fine_dx_nm) + 1) * fine_dx_nm
        self.x_coarse = np.arange(-round(coarse_half_nm / coarse_dx_nm),
                                  round(coarse_half_nm / coarse_dx_nm) + 1) * coarse_dx_nm
        self.dx_fine, self.dx_coarse = fine_dx_nm, coarse_dx_nm
        self.K_fine = np.exp(2j * np.pi * np.outer(self.x_fine, f)).astype(np.complex64)
        self.K_coarse = np.exp(2j * np.pi * np.outer(self.x_coarse, f)).astype(np.complex64)
        self.K_all = np.vstack([self.K_fine, self.K_coarse])
        self.fx = (2j * np.pi * f).astype(np.complex64)[None, :]  # d/dx kernel factor (columns = x)
        self._pupil_power = float(np.sum((self.amplitude * self.aperture) ** 2)) * self.df**2

        # metric helpers
        n = self.x_fine.size
        c = n // 2
        self._sr, self._sc, self._w, self._r2 = _coverage_window((n, n), fine_dx_nm, (c, c), r_opt_nm)
        Xf, Yf = np.meshgrid(self.x_fine, self.x_fine, indexing="xy")
        Xc2, Yc2 = np.meshgrid(self.x_coarse, self.x_coarse, indexing="xy")
        self._rc = np.hypot(Xc2, Yc2)
        self._rf = np.hypot(Xf, Yf)
        self._coarse_disk = self._rc <= total_radius_nm
        self._order_f = np.argsort(self._rf.ravel())
        self._order_c = np.argsort(self._rc.ravel())
        self.photometry = photometry or Photometry(camera_pixel_nm=sysm.camera_pixel_object_nm)
        p = self.photometry
        self._W = _overlap_matrix(n, fine_dx_nm, p.roi_px, p.camera_pixel_nm)
        self.reference_peak = None
        self.reference_power = None
        ref = self.plane(0.0, np.zeros_like(self.amplitude), ideal_reference=True)
        self.reference_peak = ref["peak"]
        self.reference_power = ref["power_ropt"]

    # ------------------------------------------------------------------ #
    def displayed_phase(self, mask: str | None, params: Mapping | None = None) -> np.ndarray:
        """Phase (waves) seen by the light on the sub-sample grid."""
        d = self.device
        if mask is None or mask == "none":
            return np.zeros(self.aperture.shape, dtype=np.float32)
        if d.ideal:
            ph = build_mask(mask, self.mask_grid, self.system, params)
            return ph.astype(np.float32)
        ph = build_mask(mask, self.pixel_grid, self.system, params)  # value at pixel centres
        if d.wrap:
            ph = np.mod(ph, 1.0)
        if d.bit_depth:
            L = 2 ** int(d.bit_depth)
            ph = np.mod(np.round(ph * L), L) / L if d.wrap else np.round(ph * L) / L
        ph = np.kron(ph, np.ones((self.s, self.s)))  # each pixel -> s x s samples
        if d.fringe_sigma_px > 0:
            ph = gaussian_filter(ph, sigma=d.fringe_sigma_px * self.s, mode="nearest")
        return ph.astype(np.float32)

    def slm_levels(self, mask: str | None, params: Mapping | None = None) -> np.ndarray:
        """The 8-bit image to send to the SLM (2 R_mask square), for inspection."""
        if mask is None or mask == "none":
            return np.zeros((2 * self.half_px,) * 2, dtype=np.uint16)
        ph = np.mod(build_mask(mask, self.pixel_grid, self.system, params), 1.0)
        L = 2 ** int(self.device.bit_depth or 8)
        return (np.mod(np.round(ph * L), L)).astype(np.uint16)

    # ------------------------------------------------------------------ #
    def _fields(self, z_nm: float, phase: np.ndarray, deriv: bool):
        P = self.amplitude * np.exp(2j * np.pi * (phase + np.float32(z_nm) * self.defocus_per_nm))
        P = np.where(self.aperture, P, 0).astype(np.complex64)
        T = self.K_all @ P  # rows: y (fine then coarse), cols: pupil x
        nf = self.x_fine.size
        Ef = T[:nf] @ self.K_fine.T
        Ec = T[nf:] @ self.K_coarse.T
        dEf = (T[:nf] * self.fx) @ self.K_fine.T if deriv else None
        return Ef, Ec, dEf

    def plane(self, z_nm: float, phase: np.ndarray, with_crlb: bool = True,
              keep_images: bool = False, ideal_reference: bool = False) -> dict:
        Ef, Ec, dEf = self._fields(z_nm, phase, deriv=with_crlb and not ideal_reference)
        norm = self.df**4 / self._pupil_power
        If = (np.abs(Ef) ** 2).astype(np.float64) * norm  # 1/nm^2, unit total power
        Ic = (np.abs(Ec) ** 2).astype(np.float64) * norm
        a_f, a_c = self.dx_fine**2, self.dx_coarse**2
        patch = If[self._sr, self._sc] * a_f
        p_opt = float((patch * self._w).sum())
        m2 = float((patch * self._w * self._r2).sum() / p_opt) if p_opt > 0 else np.nan
        p_tot = float(Ic[self._coarse_disk].sum() * a_c)
        out = {"z_nm": float(z_nm), "peak": float(If.max()), "power_ropt": p_opt, "total_power": p_tot,
               "encircled_energy": p_opt / p_tot, "second_moment_nm2": m2}
        if ideal_reference:
            return out
        out["strehl"] = out["peak"] / self.reference_peak
        c = self.x_fine.size // 2
        out["strehl_axial"] = float(If[c, c]) / self.reference_peak
        out["power_weighted_m2"] = m2 * self.reference_power / p_opt
        # EE50: fine grid if it holds half the light, else coarse
        cf = np.cumsum(If.ravel()[self._order_f]) * a_f
        if cf[-1] >= 0.5:
            j = int(np.searchsorted(cf, 0.5))
            out["ee50_nm"] = float(self._rf.ravel()[self._order_f][j])
        else:
            cc = np.cumsum(Ic.ravel()[self._order_c]) * a_c
            j = int(np.searchsorted(cc, 0.5))
            out["ee50_nm"] = float(self._rc.ravel()[self._order_c][min(j, cc.size - 1)])
        if with_crlb:
            ph = self.photometry
            dIf = 2.0 * np.real(np.conj(Ef) * dEf).astype(np.float64) * norm
            mu_s = self._W @ (If * a_f) @ self._W.T
            dmu_s = self._W @ (dIf * a_f) @ self._W.T
            out["crlb_nm"] = crlb_from(mu_s, dmu_s, ph)
            out["crlb_highbg_nm"] = crlb_from(mu_s, dmu_s, replace(ph, background=50.0))
        if keep_images:
            out["image_fine"] = If.astype(np.float32)
            out["image_camera"] = (self._W @ (If * a_f) @ self._W.T).astype(np.float32)
            out["row_coarse"] = Ic[Ic.shape[0] // 2].astype(np.float32)
        return out

    def sweep(self, z_nm: Sequence[float], phase: np.ndarray, with_crlb: bool = True,
              keep_rows: bool = True) -> dict:
        rows = [self.plane(z, phase, with_crlb=with_crlb, keep_images=keep_rows) for z in z_nm]
        keys = [k for k in rows[0] if not k.startswith("image") and not k.startswith("row")]
        out = {k: np.array([r[k] for r in rows]) for k in keys}
        if keep_rows:
            nf = self.x_fine.size
            out["xz_fine"] = np.array([r["image_fine"][nf // 2] for r in rows])
            out["xz_coarse"] = np.array([r["row_coarse"] for r in rows])
        return out


def crlb_from(mu_signal: np.ndarray, dmu_signal: np.ndarray, ph: Photometry) -> float:
    """Lateral CRLB (nm) from per-pixel expected signal fraction and its x-derivative."""
    mu = ph.n_photons * mu_signal + ph.background
    dmu = ph.n_photons * dmu_signal
    fisher = float(np.sum(dmu**2 / (mu + ph.read_noise_e**2)))
    return 1.0 / np.sqrt(fisher) if fisher > 0 else np.inf


def range_summary_slm(sw: dict, z_range=(-3000.0, 3000.0), s_min: float = 0.03, p: float = -10.0,
                      budget_nm: float | None = None) -> dict:
    z = sw["z_nm"]
    sel = (z >= z_range[0] - 1e-9) & (z <= z_range[1] + 1e-9)
    zz, ee, st = z[sel], sw["encircled_energy"][sel], sw["strehl"][sel]
    w = trapezoid_weights(zz)
    core = sel & (np.abs(z) <= 2000 + 1e-9)
    out = {
        "ee_softmin": soft_min(ee, p, w), "ee_min": float(ee.min()), "ee_mean": float(np.sum(w * ee) / w.sum()),
        "ee_flatness": float(ee.min() / (np.sum(w * ee) / w.sum())),
        "strehl_min": float(st.min()), "strehl_mean": float(np.sum(w * st) / w.sum()),
        "feasible": bool(st.min() >= s_min),
        "ee50_max_core_nm": float(sw["ee50_nm"][core].max()),
        "pwm2_mean_nm2": float(np.mean(sw["power_weighted_m2"][sel])),
    }
    if budget_nm:
        out["eta"] = float(np.trapezoid(sw["power_ropt"][sel], zz) / budget_nm)
    if "crlb_nm" in sw:
        c, ch = sw["crlb_nm"][sel], sw["crlb_highbg_nm"][sel]
        out.update({"crlb_max_nm": float(c.max()), "crlb_mean_nm": float(c.mean()),
                    "crlb_min_nm": float(c.min()),
                    "z_crlb_worst_nm": float(zz[np.argmax(c)]), "z_crlb_best_nm": float(zz[np.argmin(c)]),
                    "crlb_max_highbg_nm": float(ch.max()), "crlb_mean_highbg_nm": float(ch.mean())})
    return out
