"""Mask export for the microscope (psfsim/export.py, scripts/export_masks.py)."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from psfsim import OpticalSystem, build_mask
from psfsim.export import (FORMAT, export_mask, load_mask_file, parametric_section,
                           raster_grid, rasterize, to_psfsim)
from psfsim.families import (AnnularBinary, LayerCake, SphericalDefocus, ZernikeM0, ZonedFocus,
                             ZonedZernike2, family_by_name)

ROOT = Path(__file__).resolve().parents[1]
SYSTEM = OpticalSystem()
ZONED = {"radii": (0.48, 0.81), "depths_nm": (165.0, -1592.0, 2351.0),
         "pistons_waves": (0.0, 1.75, 0.84)}


def test_raster_grid_spans_the_pupil():
    g = raster_grid(8)
    np.testing.assert_allclose(g.x[0], (2 * np.arange(8) + 1) / 8 - 1)   # x along columns
    np.testing.assert_allclose(g.y[:, 0], g.x[0])                        # y along rows
    assert g.mask.sum() == np.sum(g.x**2 + g.y**2 <= 1)


@pytest.mark.parametrize("mask,params", [
    ("zoned_focus", ZONED),
    ("zoned_focus", {**ZONED, "sweeps_nm": (2000.0, -500.0, 0.0)}),
    ("zernike", {"coefficients": {4: 0.37, 11: 0.29, (3, -1): 0.1}}),
    ("spherical", {"strength_waves": 0.3}),
    ("annular", {"radii": (0.52, 0.66, 0.79), "phases_waves": (0.0, 0.5, 0.0, 0.5)}),
])
def test_parametric_roundtrip_is_exact(tmp_path, mask, params):
    """The exported parameters rebuild the same phase, and the raster is that phase."""
    path = export_mask(tmp_path, "m", mask, params, SYSTEM, raster_pixels=96)
    doc = load_mask_file(path)
    assert doc["format"] == FORMAT and doc["name"] == "m"
    assert doc["optical_system"] == {"na": 1.1, "n_medium": 1.33, "wavelength_nm": 525.0}
    grid = raster_grid(96)
    want = build_mask(mask, grid, SYSTEM, params)
    again, again_params = to_psfsim(doc["parametric"])
    np.testing.assert_allclose(build_mask(again, grid, SYSTEM, again_params), want, atol=1e-12)
    assert doc["raster_waves"].shape == (96, 96) and doc["raster_waves"].dtype == np.float32
    np.testing.assert_allclose(doc["raster_waves"], want, atol=1e-5)
    assert np.all(doc["raster_waves"][~grid.mask] == 0)


def test_raster_is_unwrapped():
    r = rasterize("zoned_focus", ZONED, SYSTEM, 64)
    assert r.max() - r.min() > 1.0  # more than one wave: not wrapped


def test_mask_without_parametric_form_is_raster_only(tmp_path):
    assert parametric_section("cubic", {"alpha_waves": 2.0}) is None
    doc = load_mask_file(export_mask(tmp_path, "cubic2", "cubic", {"alpha_waves": 2.0}, raster_pixels=32))
    assert "parametric" not in doc and doc["raster"] == {"file": "cubic2.npy", "pixels": 32}
    assert doc["source"]["psfsim_params"] == {"alpha_waves": 2.0}
    assert doc["raster_waves"][0, -1] == 0 and abs(doc["raster_waves"][16, 28]) > 0.5  # x^3 + y^3
    with pytest.raises(ValueError, match="needs a raster"):
        export_mask(tmp_path, "cubic3", "cubic", raster_pixels=None)


def test_no_raster_and_bad_names(tmp_path):
    path = export_mask(tmp_path, "z3", "zoned_focus", ZONED, raster_pixels=None, label="Three zones")
    doc = json.loads(path.read_text())
    assert "raster" not in doc and not (tmp_path / "z3.npy").exists() and doc["label"] == "Three zones"
    assert doc["parametric"]["sweeps_nm"] == [0.0, 0.0, 0.0]  # defaults written out
    with pytest.raises(ValueError, match="letters, digits and underscores"):
        export_mask(tmp_path, "bad name+", "zoned_focus", ZONED)


@pytest.mark.parametrize("fam", [SphericalDefocus(), ZernikeM0(), ZonedZernike2(), AnnularBinary(3),
                                 LayerCake(5), ZonedFocus(3), ZonedFocus(4, sweep=True)])
def test_family_by_name(fam):
    again = family_by_name(fam.name)
    assert type(again) is type(fam) and again.names == fam.names
    with pytest.raises(KeyError):
        family_by_name("annular_binary_3_sweep")


def test_export_script_writes_the_candidates(tmp_path):
    """The three bench candidates, straight from docs/phase2/results.json, into an
    exports folder and a controller checkout's mask library."""
    controller = tmp_path / "controller"
    (controller / "config").mkdir(parents=True)
    out = tmp_path / "exports"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "export_masks.py"), "--out", str(out),
                    "--raster-pixels", "64", "--controller", str(controller)], check=True)
    names = {"zoned_focus_3", "zoned_focus_4_sweep", "spherical_defocus"}
    assert {p.stem for p in out.glob("*.json")} == names == {p.stem for p in out.glob("*.npy")}
    library = controller / "config" / "mask_library"
    assert {p.name for p in library.iterdir()} == {f"{n}.json" for n in names}  # no rasters

    results = json.loads((ROOT / "docs" / "phase2" / "results.json").read_text())["configs"]
    doc = load_mask_file(out / "zoned_focus_4_sweep.json")
    x = results["zoned_focus_4_sweep"]["x"]
    assert doc["parametric"]["depths_nm"] == x[4:8] and doc["parametric"]["sweeps_nm"] == x[11:15]
    assert doc["source"]["metrics"]["ee_softmin"] == results["zoned_focus_4_sweep"]["ee_softmin"]
    doc = load_mask_file(out / "spherical_defocus.json")
    assert doc["parametric"] == {"type": "zernike", "modes": [
        {"n": 2, "m": 0, "waves_rms": results["spherical+defocus"]["x"][0]},
        {"n": 4, "m": 0, "waves_rms": results["spherical+defocus"]["x"][1]}]}


def test_committed_exports_match_the_results():
    """exports/*.json are the masks docs/phase2/results.json reports."""
    results = json.loads((ROOT / "docs" / "phase2" / "results.json").read_text())["configs"]
    for name, config in (("zoned_focus_3", "zoned_focus_3"), ("zoned_focus_4_sweep", "zoned_focus_4_sweep"),
                         ("spherical_defocus", "spherical+defocus")):
        doc = json.loads((ROOT / "exports" / f"{name}.json").read_text())
        mask, params = family_by_name(config).mask(results[config]["x"])
        assert doc["parametric"] == parametric_section(mask, params)
