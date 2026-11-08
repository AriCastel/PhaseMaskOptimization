"""Phase 1: objective, soft-min and CRLB checks."""

import numpy as np
import pytest

from psfsim import OpticalSystem, PSFSimulator, RadialPSF
from psfsim.crlb import LateralCRLB, Photometry
from psfsim.detector import _overlap_matrix
from psfsim.objective import Evaluator, Targets, soft_min


def test_soft_min_bounds():
    v = np.array([0.1, 0.2, 0.25, 0.3])
    sm = soft_min(v, -10)
    assert v.min() < sm < 1 / np.mean(1 / v)
    assert soft_min(v, -200) == pytest.approx(v.min(), rel=0.01)
    assert soft_min(np.full(5, 0.2)) == pytest.approx(0.2)


@pytest.fixture(scope="module")
def ev():
    return Evaluator(RadialPSF(), Targets())


def test_evaluator_consistent_with_phase0(ev):
    m = ev.evaluate_mask("spherical", {"strength_waves": 1 / np.pi})
    assert m["eta"] == pytest.approx(0.616, abs=0.005)  # Phase 0 table
    assert m["ee_min"] <= m["ee_softmin"] <= m["ee_mean"]
    assert 0 < m["ee_flatness"] <= 1


def test_constraints(ev):
    t = Targets(s_min=0.0, r_max_nm=np.inf)
    e2 = Evaluator(ev.rad, t)
    m = e2.evaluate_mask(None)
    assert m["feasible"] and m["score"] == pytest.approx(m["ee_softmin"])
    m = ev.evaluate_mask(None)  # default targets: no mask breaks the Strehl floor (S_min 0.03)
    assert not m["feasible"] and m["score"] < m["ee_softmin"]


@pytest.mark.parametrize("mask,z", [((None, None), 0.0), (("spherical", {"strength_waves": 1 / np.pi}), 1500.0)])
def test_crlb_matches_2d_bruteforce(mask, z):
    """Radial-model CRLB vs Fisher information of the shifted, camera-binned 2D PSF."""
    s = OpticalSystem(slm_bit_depth=None, grid_size=2048)
    sim, rad = PSFSimulator(s), RadialPSF(s)
    ph = Photometry(camera_pixel_nm=s.camera_pixel_object_nm)
    n, dx = s.grid_size, s.sim_pixel_nm
    w = _overlap_matrix(n, dx, ph.roi_px, ph.camera_pixel_nm)
    psf = sim.simulate(z_nm=z, mask=mask[0], mask_params=mask[1], slm_pixelation=False).psf
    f = np.fft.fftfreq(n, dx)
    F = np.fft.fft2(np.fft.ifftshift(psf))

    def cam_shift(x0):
        p = np.fft.fftshift(np.real(np.fft.ifft2(F * np.exp(-2j * np.pi * f[None, :] * x0))))
        return w @ p @ w.T

    mu = ph.n_photons * cam_shift(0.0) + ph.background
    dmu = ph.n_photons * (cam_shift(2.0) - cam_shift(-2.0)) / 4.0
    sig_2d = 1 / np.sqrt(np.sum(dmu**2 / (mu + ph.read_noise_e**2)))
    sig_r = LateralCRLB(ph).sigma_x(rad.r_nm, rad.simulate([z], *mask).intensity)[0]
    assert sig_r == pytest.approx(sig_2d, rel=0.02)
