"""
The forward pipeline.

    pupil amplitude  A(rho)                          [pupil.py]
        x exp(i 2pi [ W_aberration + W_defocus(z) + W_SLM ])
                                                     [zernike.py, pupil.py, masks.py, slm.py]
    -> P, the generalised pupil function
    -> PSF = |FT{P}|^2                               (fine grid)
    -> camera image = area-integrate(PSF)            [detector.py]

Two branches of phase are kept separate on purpose:

* the **system** branch (aberrations + depth-induced defocus) is continuous —
  it is what the sample and optics do to the wavefront;
* the **SLM** branch (engineered mask + AO correction) is discretised, because
  that is the only phase you actually control, and you control it through a
  pixelated, quantised device.

Keeping them apart is what makes the AO problem well-posed later: the loop
searches over the SLM branch to cancel the system branch, and the simulation
honestly reports what the device can and cannot represent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping, Optional, Sequence

import numpy as np

from .config import OpticalSystem
from .detector import CameraImage, downsample_to_camera
from .masks import build_mask
from .pupil import (
    PupilGrid,
    defocus_phase,
    make_pupil_grid,
    peak_to_valley,
    pupil_amplitude,
    rms_over_pupil,
)
from .slm import SLMResult, apply_slm
from .zernike import ModeSpec, ZernikeBasis, noll_range, resolve_mode


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class PSFResult:
    """
    Every intermediate of one simulation run.

    Phase maps are in **waves** (multiply by 2*pi for radians). The PSF is
    normalised to unit total energy on the fine grid, so the camera image is in
    the same units and `strehl` is directly comparable between runs.
    """

    system: OpticalSystem
    grid: PupilGrid

    amplitude: np.ndarray
    phase_aberration: np.ndarray
    phase_defocus: np.ndarray
    phase_mask_ideal: np.ndarray
    phase_mask_applied: np.ndarray
    pupil: np.ndarray  # complex generalised pupil function

    psf: np.ndarray  # fine grid, unit total energy
    camera: CameraImage
    slm: SLMResult

    z_nm: float = 0.0
    aberrations: Dict[int, float] = field(default_factory=dict)
    mask_name: Optional[str] = None
    _reference_peak: Optional[float] = None

    # -- phase summaries ----------------------------------------------------
    @property
    def phase_system(self) -> np.ndarray:
        """Aberration + depth defocus, i.e. what AO would have to cancel."""
        return self.phase_aberration + self.phase_defocus

    @property
    def phase_total(self) -> np.ndarray:
        """The phase actually imposed on the pupil, including the SLM's 2*pi
        wrapping. This is what forms the PSF."""
        return np.where(
            self.grid.mask, self.phase_system + self.phase_mask_applied, 0.0
        )

    @property
    def phase_total_unwrapped(self) -> np.ndarray:
        """
        The equivalent continuous wavefront: system phase plus the *ideal*
        (unwrapped, unpixelated) mask.

        Wavefront statistics must be computed on this, not on `phase_total`: the
        RMS of a wrapped map is dominated by the 2*pi seams and is meaningless —
        a perfectly corrected wavefront displayed through a wrapping SLM would
        report ~0.3 waves RMS while giving Strehl 1.
        """
        return np.where(
            self.grid.mask, self.phase_system + self.phase_mask_ideal, 0.0
        )

    @property
    def wavefront_rms_waves(self) -> float:
        """RMS of the equivalent continuous wavefront, piston removed."""
        return rms_over_pupil(self.phase_total_unwrapped, self.grid)

    @property
    def wavefront_pv_waves(self) -> float:
        return peak_to_valley(self.phase_total_unwrapped, self.grid)

    # -- image-plane coordinates -------------------------------------------
    @property
    def x_nm(self) -> np.ndarray:
        n = self.system.grid_size
        return (np.arange(n) - n // 2) * self.system.sim_pixel_nm

    @property
    def extent_nm(self) -> tuple[float, float, float, float]:
        """(left, right, bottom, top) for imshow of the fine PSF."""
        h = self.system.field_of_view_nm / 2
        return (-h, h, -h, h)

    # -- transfer functions -------------------------------------------------
    @property
    def otf(self) -> np.ndarray:
        """
        Incoherent OTF, normalised to OTF(0) = 1, fftshifted (DC at centre).

        Derived from the PSF rather than assembled in the pupil plane, because
        the OTF is the *autocorrelation* of the pupil function, not the pupil
        function itself — aberrations and masks are not additive on it.
        """
        otf = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(self.psf)))
        return otf / np.abs(otf).max()

    @property
    def mtf(self) -> np.ndarray:
        return np.abs(self.otf)

    @property
    def ptf(self) -> np.ndarray:
        """Phase transfer function, in radians. Non-zero PTF means the mask or
        aberration shifts image features, which matters for localisation bias."""
        return np.angle(self.otf)

    @property
    def frequencies_cycles_per_um(self) -> np.ndarray:
        n = self.system.grid_size
        df = 1.0 / (n * self.system.sim_pixel_nm)
        return (np.arange(n) - n // 2) * df * 1000.0

    def radial_mtf(self, n_bins: int = 200) -> tuple[np.ndarray, np.ndarray]:
        """Azimuthally averaged MTF vs frequency in cycles/um."""
        f = self.frequencies_cycles_per_um
        fx, fy = np.meshgrid(f, f, indexing="xy")
        fr = np.hypot(fx, fy)
        cutoff = self.system.cutoff_frequency * 1000.0
        edges = np.linspace(0, cutoff, n_bins + 1)
        idx = np.clip(np.digitize(fr.ravel(), edges) - 1, 0, n_bins - 1)
        keep = fr.ravel() <= cutoff
        sums = np.bincount(idx[keep], weights=self.mtf.ravel()[keep], minlength=n_bins)
        counts = np.bincount(idx[keep], minlength=n_bins)
        with np.errstate(invalid="ignore"):
            profile = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)
        return 0.5 * (edges[:-1] + edges[1:]), profile

    # -- quality metrics ----------------------------------------------------
    @property
    def strehl(self) -> float:
        """Peak intensity relative to the unaberrated, unmasked PSF."""
        if self._reference_peak is None:
            return float("nan")
        return float(self.psf.max() / self._reference_peak)

    @property
    def strehl_marechal(self) -> float:
        """exp(-(2*pi*sigma)^2) estimate. Only valid for small aberrations —
        comparing it against `strehl` is a quick sanity check on whether you are
        still in the weakly aberrated regime (they diverge above ~0.2 waves)."""
        return float(np.exp(-((2 * np.pi * self.wavefront_rms_waves) ** 2)))

    @property
    def encircled_energy_radius_nm(self) -> float:
        """Radius containing 50% of the fine-grid energy — a mask-agnostic size
        metric that, unlike FWHM, behaves sensibly for engineered PSFs."""
        n = self.system.grid_size
        x = self.x_nm
        r = np.hypot(*np.meshgrid(x, x, indexing="xy"))
        order = np.argsort(r.ravel())
        cum = np.cumsum(self.psf.ravel()[order])
        return float(r.ravel()[order][np.searchsorted(cum, 0.5 * cum[-1])])

    @property
    def edge_energy_fraction(self) -> float:
        """
        Fraction of energy in the outer 5% border of the fine grid. If this is
        not small, the PSF is wrapping around the FFT window and the result is
        aliased — increase `pupil_radius_px` (which increases the FOV).
        """
        n = self.system.grid_size
        b = max(1, int(0.05 * n))
        interior = self.psf[b:-b, b:-b].sum()
        return float(1.0 - interior / self.psf.sum())

    def diagnostics(self) -> str:
        lines = [
            f"z = {self.z_nm:+.0f} nm, mask = {self.mask_name or 'none'}",
            f"  wavefront RMS (total) : {self.wavefront_rms_waves:.3f} waves",
            f"  wavefront PV (total)  : {self.wavefront_pv_waves:.3f} waves",
            f"  Strehl (from PSF)     : {self.strehl:.4f}",
            f"  50% energy radius     : {self.encircled_energy_radius_nm:.0f} nm",
            f"  camera image          : {self.camera.shape[0]}^2 px "
            f"@ {self.camera.pixel_size_nm:.1f} nm",
            f"  SLM samples/pixel     : {self.slm.samples_per_slm_pixel:.2f}"
            f"{'' if self.slm.resolves_device_pixels else '  [device finer than grid]'}"
            f"{'  [MASK UNDERSAMPLED]' if self.slm.undersampled else ''}",
            f"  edge energy           : {self.edge_energy_fraction*100:.2f} %"
            f"{'  [ALIASING]' if self.edge_energy_fraction > 0.01 else ''}",
        ]
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Simulator
# --------------------------------------------------------------------------- #
class PSFSimulator:
    """
    Forward model for the microscope PSF.

    The grid, the pupil amplitude, the Zernike basis, and the diffraction-limited
    reference peak are all built once and reused, so sweeping depth or searching
    over coefficients is cheap.

    Examples
    --------
    >>> sim = PSFSimulator(OpticalSystem(na=1.0, wavelength_nm=561))
    >>> r = sim.simulate(z_nm=500, aberrations={11: 0.2},
    ...                  mask="spherical", mask_params={"strength_waves": 1.5})
    >>> r.camera.image.shape, r.strehl
    """

    def __init__(
        self,
        system: OpticalSystem | None = None,
        basis_modes: Iterable[ModeSpec] | None = None,
    ) -> None:
        self.system = system or OpticalSystem()
        self.grid = make_pupil_grid(self.system)
        self.amplitude = pupil_amplitude(self.grid, self.system)
        self.basis = ZernikeBasis(
            basis_modes if basis_modes is not None else noll_range(21),
            self.grid.rho,
            self.grid.theta,
            self.grid.mask,
        )
        self._reference_peak = self._psf_from_pupil(
            self.amplitude.astype(complex)
        ).max()

    # -- core ---------------------------------------------------------------
    def _psf_from_pupil(self, pupil: np.ndarray) -> np.ndarray:
        field = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(pupil)))
        psf = np.abs(field) ** 2
        return psf / psf.sum()

    def wavefront_from_coefficients(
        self, coefficients: Mapping[ModeSpec, float] | None
    ) -> np.ndarray:
        """Zernike coefficients (waves RMS) -> 2D phase map, via the cached basis.

        Modes outside the cached basis are evaluated on the fly, so you are not
        limited to the first 21 Noll terms.
        """
        if not coefficients:
            return np.zeros_like(self.grid.rho)
        from .zernike import zernike_map

        out = np.zeros_like(self.grid.rho)
        for mode, coeff in coefficients.items():
            nm = resolve_mode(mode)
            if nm in self.basis.modes:
                out += coeff * self.basis[nm]
            else:
                out += coeff * zernike_map(
                    nm, self.grid.rho, self.grid.theta, self.grid.mask
                )
        return out

    def simulate(
        self,
        z_nm: float = 0.0,
        aberrations: Mapping[ModeSpec, float] | None = None,
        mask: str | None = None,
        mask_params: Mapping | None = None,
        correction: Mapping[ModeSpec, float] | None = None,
        extra_pupil_phase: np.ndarray | None = None,
        exact_defocus: bool = True,
        slm_pixelation: bool = True,
        slm_sampling: str = "center",
        camera_pixel_nm: float | None = None,
        camera_fov_nm: float | None = None,
    ) -> PSFResult:
        """
        Run one full forward pass.

        Parameters
        ----------
        z_nm : axial position of the point source relative to the focal plane.
        aberrations : {Noll index or (n, m): RMS waves} — the *system* aberration
            you are trying to overcome.
        mask, mask_params : engineered phase mask on the SLM, e.g.
            mask="spherical", mask_params={"strength_waves": 1.5}.
        correction : Zernike coefficients added **on the SLM**, alongside the
            mask. This is the AO hook: an AO loop optimises this dict, and
            because it goes through the device model it is automatically subject
            to pixelation and quantisation, as it would be on the bench.
        extra_pupil_phase : arbitrary phase map in waves, added to the system
            branch. Use for measured wavefronts or non-Zernike sample-induced
            phase (e.g. a refractive-index boundary).
        camera_pixel_nm : override the object-side camera pitch for this run.
        camera_fov_nm : crop the camera image to this field of view (object-side
            nm) around the axis. The full simulated field is much larger than a
            useful ROI, so cropping keeps the main output a manageable size.
        """
        aberr = self.wavefront_from_coefficients(aberrations)
        if extra_pupil_phase is not None:
            aberr = aberr + np.where(self.grid.mask, extra_pupil_phase, 0.0)

        defocus = defocus_phase(self.grid, self.system, z_nm, exact=exact_defocus)

        mask_ideal = build_mask(mask, self.grid, self.system, mask_params)
        if correction:
            mask_ideal = mask_ideal + self.wavefront_from_coefficients(correction)

        slm = apply_slm(
            mask_ideal,
            self.grid,
            self.system,
            mode=slm_sampling,
            pixelate_enabled=slm_pixelation,
        )

        total = np.where(self.grid.mask, aberr + defocus + slm.applied_waves, 0.0)
        pupil = self.amplitude * np.exp(2j * np.pi * total)
        psf = self._psf_from_pupil(pupil)

        camera = downsample_to_camera(
            psf,
            self.system.sim_pixel_nm,
            camera_pixel_nm or self.system.camera_pixel_object_nm,
            fov_nm=camera_fov_nm,
        )

        return PSFResult(
            system=self.system,
            grid=self.grid,
            amplitude=self.amplitude,
            phase_aberration=aberr,
            phase_defocus=defocus,
            phase_mask_ideal=mask_ideal,
            phase_mask_applied=slm.applied_waves,
            pupil=pupil,
            psf=psf,
            camera=camera,
            slm=slm,
            z_nm=z_nm,
            aberrations=dict(aberrations or {}),
            mask_name=mask,
            _reference_peak=self._reference_peak,
        )

    # -- convenience --------------------------------------------------------
    def iter_stack(self, z_values: Sequence[float], **kwargs):
        """
        Generator over depth. Each `PSFResult` holds several full-grid arrays
        (~60 MB at a 1024^2 grid), so a depth sweep that keeps every result
        alive will exhaust memory. Iterate and extract what you need instead;
        `simulate_stack` is only for short sweeps you really want to retain.
        """
        kwargs.pop("z_nm", None)
        for z in z_values:
            yield self.simulate(z_nm=float(z), **kwargs)

    def simulate_stack(self, z_values: Sequence[float], **kwargs) -> list[PSFResult]:
        """Same settings, swept over depth, all results retained. Keep the
        number of planes modest — see `iter_stack` for the memory caveat."""
        return list(self.iter_stack(z_values, **kwargs))

    def axial_profile(
        self, z_values: Sequence[float], **kwargs
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        (z, Strehl) over a depth range — the one-number-per-plane summary used to
        compare EDOF masks. A good mask trades a lower peak Strehl for a much
        flatter curve, so what matters is the width of this curve relative to
        `system.depth_of_field_nm`, not its height.
        """
        strehl = [r.strehl for r in self.iter_stack(z_values, **kwargs)]
        return np.asarray(z_values, dtype=float), np.array(strehl)

    def xz_slice(self, z_values: Sequence[float], **kwargs) -> np.ndarray:
        """(n_z, n_x) axial section through the fine PSF, for x-z display."""
        row = self.system.grid_size // 2
        return np.array([r.psf[row, :].copy() for r in self.iter_stack(z_values, **kwargs)])

    def camera_stack(
        self, z_values: Sequence[float], camera_fov_nm: float = 4000.0, **kwargs
    ) -> np.ndarray:
        """
        (n_z, n_y, n_x) stack of camera images — the in-silico z-stack that an
        AO metric or a localisation algorithm would be fed.
        """
        kwargs.pop("camera_fov_nm", None)
        return np.array(
            [r.camera.image.copy() for r in self.iter_stack(z_values, camera_fov_nm=camera_fov_nm, **kwargs)]
        )

    def depth_metrics(self, z_values: Sequence[float], **kwargs) -> dict:
        """
        Per-plane scalar metrics over a depth sweep, without retaining the full
        results: Strehl (peak relative to the diffraction limit), 50%
        encircled-energy radius in nm, and total wavefront RMS in waves.
        Same kwargs as `simulate`.
        """
        out = {"strehl": [], "ee50_nm": [], "rms_waves": []}
        for r in self.iter_stack(z_values, **kwargs):
            out["strehl"].append(r.strehl)
            out["ee50_nm"].append(r.encircled_energy_radius_nm)
            out["rms_waves"].append(r.wavefront_rms_waves)
        return {k: np.asarray(v) for k, v in out.items()}

    def psf_volume(
        self,
        z_values: Sequence[float],
        fov_nm: float = 4000.0,
        use_camera: bool = False,
        **kwargs,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        (n_z, ny, nx) volume of same-FOV, same-sampling PSF crops across depth —
        the data behind a top-down / multi-slice PSF-engineering figure (see
        `psfsim.visualize.plot_topdown_with_slices`): each plane is a "front"
        (x, y) view at one z, and taking a row through the stack at fixed y
        gives the "top-down" (z, x) view.

        use_camera=False (default) crops the fine simulation grid, which is the
        resolution you want when judging a mask's shape rather than what the
        camera pixel grid does to it; use_camera=True crops the camera-sampled
        image instead.

        Returns (volume, x_nm, z_nm) — x_nm is the shared lateral coordinate
        axis (both dimensions are square and identical), z_nm is `z_values`
        echoed back as an array, for convenient unpacking into a plotting call.
        """
        from .detector import crop_center

        volume = None
        x_nm = None
        for i, r in enumerate(self.iter_stack(z_values, **kwargs)):
            if use_camera:
                img, px = r.camera.image, r.camera.pixel_size_nm
            else:
                img, px = r.psf, self.system.sim_pixel_nm
            cropped, x_nm = crop_center(img, px, fov_nm)
            if volume is None:
                volume = np.empty((len(z_values),) + cropped.shape, dtype=cropped.dtype)
            volume[i] = cropped
        return volume, x_nm, np.asarray(z_values, dtype=float)
