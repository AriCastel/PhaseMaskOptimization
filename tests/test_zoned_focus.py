"""Phase 2: zoned_focus mask and search families."""

import numpy as np
import pytest

from psfsim import OpticalSystem, PSFSimulator, RadialPSF
from psfsim.families import AnnularBinary, LayerCake, SphericalDefocus, ZernikeM0, ZonedFocus, ZonedZernike2
from psfsim.masks import equal_dc_radii, rho_to_t, t_to_rho
from psfsim.metrics import plane_metrics

SYSTEM = OpticalSystem(slm_bit_depth=None)


@pytest.fixture(scope="module")
def rad():
    return RadialPSF(SYSTEM)


def test_equal_dc_radii():
    r = equal_dc_radii(4, SYSTEM)
    np.testing.assert_allclose(rho_to_t(np.array(r), SYSTEM), [0.25, 0.5, 0.75], atol=1e-12)
    np.testing.assert_allclose(rho_to_t(t_to_rho(np.linspace(0, 1, 7), SYSTEM), SYSTEM), np.linspace(0, 1, 7),
                               atol=1e-12)


@pytest.mark.parametrize("zk", [-2500.0, 1500.0])
def test_single_depth_moves_focus_exactly(rad, zk):
    """All zones at the same depth = exact refocus: Strehl 1 at z = z_k."""
    st = rad.simulate([zk], "zoned_focus", {"radii": (0.4, 0.8), "depths_nm": (zk, zk, zk)})
    assert st["strehl"][0] == pytest.approx(1.0, abs=1e-6)


def test_sweep_spreads_focus(rad):
    """One zone with sweep s covers depths 0 +- s/2: Strehl at +-1 um rises
    far above the unswept value. (The line focus is not symmetric: at high NA
    the energy per unit c = cos theta scales ~c^2, so the pupil centre, which a
    positive sweep sends to +z, carries ~3x more light per unit depth.)"""
    z = np.array([-1000.0, 1000.0])
    one = rad.simulate(z, "zoned_focus", {"radii": (), "depths_nm": (0.0,)})
    sw = rad.simulate(z, "zoned_focus", {"radii": (), "depths_nm": (0.0,), "sweeps_nm": (3000.0,)})
    assert np.all(sw["strehl"] > 3 * one["strehl"])
    assert sw["strehl"][1] > sw["strehl"][0]


def test_global_piston_irrelevant(rad):
    z = np.array([-1500.0, 0.0, 2000.0])
    p = {"radii": (0.5, 0.8), "depths_nm": (-1000.0, 500.0, 2000.0)}
    a = rad.simulate(z, "zoned_focus", {**p, "pistons_waves": (0.0, 0.3, 0.7)})
    b = rad.simulate(z, "zoned_focus", {**p, "pistons_waves": (0.25, 0.55, 0.95)})
    np.testing.assert_allclose(a["encircled_energy"], b["encircled_energy"], rtol=1e-10)


def test_zoned_focus_radial_matches_2d(rad):
    p = {"radii": (0.55, 0.83), "depths_nm": (500.0, -1800.0, 2100.0), "pistons_waves": (0.0, 0.6, 0.3),
         "sweeps_nm": (0.0, 1000.0, -800.0)}
    sim = PSFSimulator(SYSTEM)
    ref = plane_metrics(sim.simulate().psf, SYSTEM.sim_pixel_nm, 400.0).power_ropt
    for z in (-2000.0, 0.0, 2500.0):
        r = sim.simulate(z_nm=z, mask="zoned_focus", mask_params=p, slm_pixelation=False)
        m = plane_metrics(r.psf, SYSTEM.sim_pixel_nm, 400.0, 10_000.0, ref)
        st = rad.simulate([z], "zoned_focus", p)
        assert st["encircled_energy"][0] == pytest.approx(m.encircled_energy, rel=0.02, abs=0.002)
        assert st["strehl"][0] == pytest.approx(r.strehl, abs=0.01)


@pytest.mark.parametrize("fam", [SphericalDefocus(), ZernikeM0(), ZonedZernike2(), AnnularBinary(2),
                                 LayerCake(3), ZonedFocus(3), ZonedFocus(4, sweep=True)])
def test_families_build_valid_masks(rad, fam):
    rng = np.random.default_rng(0)
    for _ in range(5):
        x = fam.lo + rng.uniform(0, 1, fam.dim) * (fam.hi - fam.lo)
        st = rad.simulate([0.0], *fam.mask(x))
        assert np.isfinite(st["encircled_energy"][0])
        assert isinstance(fam.describe(x), str)
