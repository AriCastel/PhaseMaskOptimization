"""Zoned radial-Zernike family (free Z(n, 0) up to order 6 per zone)."""

import json
from pathlib import Path

import numpy as np
import pytest

from psfsim import OpticalSystem, RadialPSF
from psfsim.families import ZonedFocus, ZonedZernikeRadial, noll_m0
from psfsim.objective import Evaluator, Targets
from psfsim.zernike import noll_to_nm

SYSTEM = OpticalSystem(slm_bit_depth=None)
PHASE2 = json.loads((Path(__file__).resolve().parents[1] / "docs" / "phase2" / "results.json").read_text())["configs"]


@pytest.fixture(scope="module")
def ev():
    return Evaluator(RadialPSF(SYSTEM, n_pupil=1024, dr_nm=20.0), Targets())


def test_noll_m0():
    assert [noll_m0(n) for n in (0, 2, 4, 6, 8)] == [1, 4, 11, 22, 37]
    assert all(noll_to_nm(noll_m0(n)) == (n, 0) for n in (0, 2, 4, 6, 8))


@pytest.mark.parametrize("n_zones, order, dim", [(3, 6, 14), (4, 6, 19), (4, 2, 11)])
def test_dimensions_and_mask(ev, n_zones, order, dim):
    fam = ZonedZernikeRadial(n_zones, order)
    assert fam.dim == dim == len(fam.lo) == len(fam.hi)
    x = np.random.default_rng(0).uniform(fam.lo, fam.hi)
    name, p = fam.mask(x)
    assert name == "zoned_zernike" and len(p["radii"]) == n_zones - 1
    assert 1 not in p["zone_coefficients"][0]  # zone 0 piston is the global piston
    assert np.all(np.isfinite(ev.phase(name, p)))


def test_local_and_global_zernikes_span_the_same_masks(ev):
    """rho'^2 is affine in rho^2, so a zone of local Zernikes up to order 6 is a
    cubic in rho^2: a full-pupil Zernike sum of the same order is reproduced exactly."""
    fam = ZonedZernikeRadial(3, 6)
    ph = ev.phase("zernike", {"coefficients": {4: 0.4, 11: -0.3, 22: 0.1}})
    radii = (0.5, 0.8)
    x = fam.from_phase(radii, ev.rad.rho, ph)
    got = ev.phase(*fam.mask(x))
    d = np.mod(got - ph + 0.5, 1.0) - 0.5  # pistons are wrapped mod one wave
    for a, b in zip((0.0,) + radii, radii + (1.0,)):
        sel = (ev.rad.rho >= a) & (ev.rad.rho < b)
        assert np.ptp(d[sel]) < 1e-9


@pytest.mark.parametrize("n_zones, sweep", [(3, False), (4, False), (4, True)])
def test_contains_zoned_focus(ev, n_zones, sweep):
    """Order 6 reproduces the Phase 2 zoned_focus optima (exact defocus, sweep) to < 1e-3 in score."""
    src = ZonedFocus(n_zones, sweep)
    m, p = src.mask(PHASE2[src.name]["x"])
    ph = ev.phase(m, p)
    fam = ZonedZernikeRadial(n_zones, 6)
    x = fam.from_phase(p["radii"], ev.rad.rho, ph)
    assert np.all((x >= fam.lo) & (x <= fam.hi))
    a, b = ev.evaluate(ph), ev.evaluate_mask(*fam.mask(x))
    assert b["score"] == pytest.approx(a["score"], abs=1e-3)
    assert b["strehl_min"] == pytest.approx(a["strehl_min"], abs=2e-3)
