"""
Export phase masks for the microscope (the reflex3d_controller mask library).

A mask is exported as a small JSON file, `<name>.json`, optionally with a
rasterised matrix next to it, `<name>.npy`:

    {
      "format": "reflex3d_phase_mask", "version": 1,
      "name": "zoned_focus_3", "label": "...", "description": "...",
      "phase_units": "waves",
      "optical_system": {"na": 1.1, "n_medium": 1.33, "wavelength_nm": 525.0},
      "parametric": {"type": "zoned_focus", "radii": [...], "depths_nm": [...],
                     "pistons_waves": [...], "sweeps_nm": [...]},
      "raster": {"file": "zoned_focus_3.npy", "pixels": 1024},
      "source": {...}                      # provenance, free-form
    }

Parametric form (`parametric`): the mask's own parameters, which the controller
evaluates analytically on each channel's pupil, whatever its radius and centre
on the SLM. Three types are understood by the controller:

    zoned_focus   radii, depths_nm, pistons_waves, sweeps_nm  (masks.zoned_focus;
                  uses optical_system for the cos(theta) mapping)
    zernike       modes: [{"n", "m", "waves_rms"}]  (RMS-normalised; m > 0 cosine,
                  m < 0 sine - explicit (n, m), so no Noll-convention ambiguity)
    annular       radii, phases_waves              (masks.annular)

Any other registered mask is exported without a parametric section, raster only.

Raster form (`raster`): a square float32 matrix of the UNWRAPPED phase in waves,
0 outside the pupil. The pupil is the inscribed disk: element [i, j] is the
phase at normalised pupil coordinates

    x = (2 j + 1) / N - 1,    y = (2 i + 1) / N - 1        (N = pixels)

with x along SLM columns and y along SLM rows. It is unwrapped so that the
reader can resample it (nearest neighbour) before wrapping it for its own SLM.
The raster is what makes free-form masks, which have no parameters, loadable.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .config import OpticalSystem
from .masks import build_mask
from .pupil import PupilGrid
from .zernike import resolve_mode

FORMAT = "reflex3d_phase_mask"
VERSION = 1
DEFAULT_RASTER_PIXELS = 1024


def raster_grid(pixels: int) -> PupilGrid:
    """The export raster's grid: `pixels` samples across the pupil diameter."""
    c = (2.0 * np.arange(int(pixels)) + 1.0) / int(pixels) - 1.0
    x, y = np.meshgrid(c, c, indexing="xy")
    rho = np.hypot(x, y)
    return PupilGrid(x=x, y=y, rho=rho, theta=np.arctan2(y, x), mask=rho <= 1.0,
                     pupil_radius_px=pixels / 2.0)


def rasterize(mask: str, params: Mapping | None, system: OpticalSystem,
              pixels: int = DEFAULT_RASTER_PIXELS) -> np.ndarray:
    """Unwrapped phase in waves on the export raster (float32, 0 outside the pupil)."""
    return build_mask(mask, raster_grid(pixels), system, params).astype(np.float32)


def parametric_section(mask: str, params: Mapping | None) -> dict | None:
    """The `parametric` section for a registered mask, or None if the controller
    has no analytic form for it (it is then exported as a raster only)."""
    p = dict(params or {})
    if mask == "zoned_focus":
        radii = [float(r) for r in p.get("radii", (1 / np.sqrt(2),))]
        n = len(radii) + 1
        depths = [float(z) for z in p.get("depths_nm", (-1000.0, 1000.0))]
        pistons = p.get("pistons_waves")
        sweeps = p.get("sweeps_nm")
        return {
            "type": "zoned_focus",
            "radii": radii,
            "depths_nm": depths,
            "pistons_waves": [0.0] * n if pistons is None else [float(v) for v in pistons],
            "sweeps_nm": [0.0] * n if sweeps is None else [float(v) for v in sweeps],
        }
    if mask == "zernike":
        modes = []
        for mode, c in (p.get("coefficients") or {}).items():
            n, m = resolve_mode(mode)
            modes.append({"n": n, "m": m, "waves_rms": float(c)})
        return {"type": "zernike", "modes": modes}
    if mask == "spherical" and p.get("normalized", True):
        return {"type": "zernike",
                "modes": [{"n": 4, "m": 0, "waves_rms": float(p.get("strength_waves", 1.5))}]}
    if mask == "annular":
        return {
            "type": "annular",
            "radii": [float(r) for r in p.get("radii", (0.76, 0.93))],
            "phases_waves": [float(v) for v in p.get("phases_waves", (0.0, 0.5, 0.0))],
        }
    return None


def to_psfsim(parametric: Mapping) -> tuple[str, dict]:
    """(mask, params) for `build_mask` / the simulators from a `parametric` section."""
    kind = parametric["type"]
    if kind == "zoned_focus":
        return "zoned_focus", {k: tuple(parametric[k])
                               for k in ("radii", "depths_nm", "pistons_waves", "sweeps_nm")}
    if kind == "zernike":
        return "zernike", {"coefficients": {(int(d["n"]), int(d["m"])): float(d["waves_rms"])
                                            for d in parametric["modes"]}}
    if kind == "annular":
        return "annular", {"radii": tuple(parametric["radii"]),
                           "phases_waves": tuple(parametric["phases_waves"])}
    raise ValueError(f"unknown parametric mask type {kind!r}")


def _jsonable(obj):
    if isinstance(obj, Mapping):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def export_mask(
    out_dir,
    name: str,
    mask: str,
    params: Mapping | None = None,
    system: OpticalSystem | None = None,
    label: str | None = None,
    description: str = "",
    raster_pixels: int | None = DEFAULT_RASTER_PIXELS,
    source: Mapping | None = None,
) -> Path:
    """
    Write `<out_dir>/<name>.json` (and `<name>.npy` unless raster_pixels is
    None) for a registered mask. Returns the JSON path.

    `name` becomes the file name and the mask's name in the controller, so
    keep it to letters, digits and underscores.
    """
    if not name.replace("_", "").isalnum():
        raise ValueError(f"mask name {name!r} must be letters, digits and underscores only")
    system = system or OpticalSystem()
    parametric = parametric_section(mask, params)
    if parametric is None and raster_pixels is None:
        raise ValueError(f"mask {mask!r} has no parametric form for the controller; "
                         f"it needs a raster (raster_pixels)")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    doc = {
        "format": FORMAT,
        "version": VERSION,
        "name": name,
        "label": label or name,
        "description": description,
        "phase_units": "waves",
        "optical_system": {"na": system.na, "n_medium": system.n_medium,
                           "wavelength_nm": system.wavelength_nm},
    }
    if parametric is not None:
        doc["parametric"] = parametric
    if raster_pixels is not None:
        np.save(out_dir / f"{name}.npy", rasterize(mask, params, system, raster_pixels))
        doc["raster"] = {"file": f"{name}.npy", "pixels": int(raster_pixels)}
    doc["source"] = _jsonable({"psfsim_mask": mask,
                               **({} if parametric is not None else {"psfsim_params": params or {}}),
                               **(source or {})})

    path = out_dir / f"{name}.json"
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    return path


def load_mask_file(path) -> dict:
    """Read an exported mask file back. Adds `raster_waves` (the matrix) when the
    raster file is there."""
    path = Path(path)
    with open(path) as f:
        doc = json.load(f)
    if doc.get("format") != FORMAT:
        raise ValueError(f"{path} is not a {FORMAT} file")
    raster = doc.get("raster")
    if raster and (path.parent / raster["file"]).exists():
        doc["raster_waves"] = np.load(path.parent / raster["file"])
    return doc
