"""Phase 4a: pixel-level SLM simulator."""

import numpy as np
import pytest

from psfsim import OpticalSystem, RadialPSF
from psfsim.crlb import LateralCRLB
from psfsim.slm_psf import SLMDevice, SLMPSF

MASK = ("zoned_focus", {"radii": (0.55, 0.83), "depths_nm": (200.0, -1600.0, 2300.0),
                        "pistons_waves": (0.0, 0.7, 0.8)})


@pytest.fixture(scope="module")
def ideal():
    return SLMPSF(OpticalSystem(), SLMDevice(ideal=True), subsample=2)


@pytest.mark.parametrize("z", [-2500.0, 0.0, 1500.0])
def test_ideal_matches_radial(ideal, z):
    rad = RadialPSF()
    a = ideal.plane(z, ideal.displayed_phase(*MASK))
    b = rad.simulate([z], *MASK)
    assert a["strehl"] == pytest.approx(b["strehl"][0], rel=0.01, abs=1e-3)
    assert a["encircled_energy"] == pytest.approx(b["encircled_energy"][0], rel=0.01)
    c = LateralCRLB().sigma_x(rad.r_nm, b.intensity)[0]
    assert a["crlb_nm"] == pytest.approx(c, rel=0.02)


def test_cut_alone_is_exact():
    """Without quantisation or fringing the 2 pi cut changes nothing."""
    s = OpticalSystem()
    cut = SLMPSF(s, SLMDevice(bit_depth=None, wrap=True), subsample=1)
    raw = SLMPSF(s, SLMDevice(bit_depth=None, wrap=False), subsample=1)
    for z in (-2000.0, 1000.0):
        a = cut.plane(z, cut.displayed_phase(*MASK), with_crlb=False)
        b = raw.plane(z, raw.displayed_phase(*MASK), with_crlb=False)
        assert a["encircled_energy"] == pytest.approx(b["encircled_energy"], rel=1e-4)


def test_slm_levels():
    sim = SLMPSF(OpticalSystem(), SLMDevice(), subsample=1)
    lv = sim.slm_levels(*MASK)
    assert lv.shape == (724, 724) and lv.min() >= 0 and lv.max() <= 255
    assert np.all(sim.slm_levels(None) == 0)
