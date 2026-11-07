"""Phase 0 gate tests: radial model vs 2D model, budget invariance, convergence."""

import numpy as np
import pytest

from psfsim import OpticalSystem, PSFSimulator, RadialPSF
from psfsim.metrics import plane_metrics, range_summary
from psfsim.zernike import noll_to_nm

SYSTEM = OpticalSystem(slm_bit_depth=None)
SPH = ("spherical", {"strength_waves": 1 / np.pi})  # 2 rad RMS
ZONED = ("zoned_zernike", {"zone_coefficients": ({4: 0.75}, {4: -0.75})})


@pytest.fixture(scope="module")
def rad():
    return RadialPSF(SYSTEM)


def test_noll_indices():
    assert [noll_to_nm(j) for j in (1, 2, 3, 4, 5, 6, 11, 22, 37)] == [
        (0, 0), (1, 1), (1, -1), (2, 0), (2, -2), (2, 2), (4, 0), (6, 0), (8, 0)]


def test_normalisation(rad):
    st = rad.simulate([0.0])
    assert st["strehl"][0] == pytest.approx(1.0, abs=1e-12)
    assert st["total_power"][0] == pytest.approx(1.0, abs=0.01)  # light beyond 10 um is ~0.5 %


def test_rejects_non_symmetric(rad):
    with pytest.raises(ValueError):
        rad.simulate([0.0], "cubic")
    with pytest.raises(ValueError):
        rad.simulate([0.0], aberrations={7: 0.1})


@pytest.mark.parametrize("mask", [(None, None), SPH, ZONED])
@pytest.mark.parametrize("z", [-2500.0, 0.0, 1500.0])
def test_matches_2d(rad, mask, z):
    sim = PSFSimulator(SYSTEM)
    px = SYSTEM.sim_pixel_nm
    ref = plane_metrics(sim.simulate().psf, px, 400.0).power_ropt
    r = sim.simulate(z_nm=z, mask=mask[0], mask_params=mask[1], slm_pixelation=False)
    m = plane_metrics(r.psf, px, 400.0, 10_000.0, ref)
    st = rad.simulate([z], *mask)
    # Gate: 1-2 % (absolute tolerances are relative to the in-focus scale).
    assert st["strehl"][0] == pytest.approx(r.strehl, abs=0.01)
    assert st["encircled_energy"][0] == pytest.approx(m.encircled_energy, rel=0.02, abs=0.002)
    assert st["power_weighted_m2"][0] == pytest.approx(m.power_weighted_m2, rel=0.03)
    # 2D EE50 is quantised to pixel-centre radii (~60 nm grid)
    assert st["ee50_nm"][0] == pytest.approx(r.encircled_energy_radius_nm, abs=0.6 * px)


def test_converged(rad):
    fine = RadialPSF(SYSTEM, n_pupil=8192, dr_nm=5.0)
    z = [-8000.0, 0.0, 3000.0, 8000.0]
    a, b = rad.simulate(z, *SPH), fine.simulate(z, *SPH)
    for k in ("strehl", "encircled_energy", "power_ropt", "total_power"):
        np.testing.assert_allclose(a[k], b[k], rtol=1e-3, atol=1e-6)


@pytest.mark.parametrize(
    "mask",
    [(None, None), SPH, ZONED, ("annular", None),
     ("zernike", {"coefficients": {11: 0.4, 22: -0.3, 37: 0.2}})],
)
def test_strehl_budget_invariant(rad, mask):
    """Integral of on-axis Strehl over z is mask independent and equals the
    analytic value (wide window; the residual is the 1/z^2 tail)."""
    z = np.arange(-400_000.0, 400_000.0 + 1, 20.0)
    s = rad.axial_strehl(z, rad.total_phase(*mask))
    assert np.trapezoid(s, z) == pytest.approx(rad.budgets().strehl_axial_nm, rel=2e-3)


def test_power_budget_invariant(rad):
    """Integral of P(r_opt) over z: identical across masks on the same window,
    and close to the analytic infinite-range value."""
    z = np.arange(-60_000.0, 60_000.0 + 1, 25.0)
    vals = []
    for mask in [(None, None), SPH, ZONED]:
        st = rad.simulate(z, *mask)
        vals.append(np.trapezoid(st["power_ropt"], z))
    assert np.ptp(vals) / np.mean(vals) < 1e-4
    assert np.mean(vals) == pytest.approx(rad.budgets().power_ropt_nm, rel=0.01)


def test_range_summary():
    z = np.linspace(-8000, 8000, 161)
    v = np.where(np.abs(z) <= 3000, 0.2, 0.0)
    s = range_summary(z, v)
    assert s["flatness"] == pytest.approx(1.0)
    assert s["eta"] == pytest.approx(1.0, abs=0.02)  # trapezoid at the step edges
