"""
Export optimised masks for the microscope (reflex3d_controller).

Takes the best mask of each chosen configuration in docs/phase2/results.json and
writes, per mask, `exports/<name>.json` (parameters + provenance) and
`exports/<name>.npy` (rasterised phase matrix); see psfsim/export.py for the
format. By default it exports the three candidates carried to the bench:

    zoned_focus_3         3-zone zoned_focus
    zoned_focus_4_sweep   4-zone zoned_focus with per-zone sweep
    spherical_defocus     re-centred spherical (Z11 + Z4)

Usage:
    python scripts/export_masks.py
    python scripts/export_masks.py --configs annular_binary_3 zoned_focus_4
    python scripts/export_masks.py --controller ../reflex3d_controller

--controller also copies the files into that checkout's config/mask_library/,
where the controller's 'SLM phase mask' panel lists them (the .json only,
unless --controller-rasters: the controller evaluates the parameters itself).
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem  # noqa: E402
from psfsim.export import DEFAULT_RASTER_PIXELS, export_mask  # noqa: E402
from psfsim.families import family_by_name  # noqa: E402

RESULTS = ROOT / "docs" / "phase2" / "results.json"
CANDIDATES = ["zoned_focus_3", "zoned_focus_4_sweep", "spherical+defocus"]
METRICS = ["ee_softmin", "ee_min", "ee_mean", "eta", "strehl_min", "ee50_max_core_nm",
           "crlb_max_nm", "crlb_max_highbg_nm", "feasible"]

#: results.json configuration -> (file name, label shown in the controller)
NAMES = {
    "spherical+defocus": ("spherical_defocus", "Spherical + defocus (optimised)"),
    "zernike_m0": ("zernike_m0", "Zernike Z4/Z11/Z22/Z37 (optimised)"),
    "zoned_zernike2": ("zoned_zernike2", "Zoned Zernike, 2 zones"),
}


def name_and_label(config: str) -> tuple[str, str]:
    if config in NAMES:
        return NAMES[config]
    kind, n, *sweep = config.rsplit("_", 2) if config.endswith("_sweep") else (*config.rsplit("_", 1),)
    label = {"zoned_focus": f"Zoned focus, {n} zones", "annular_binary": f"Annular binary, {n} edges",
             "layer_cake": f"Layer cake, {n} zones"}.get(kind, config)
    return config, label + (" + sweep" if sweep else "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--configs", nargs="+", default=CANDIDATES,
                    help="configurations in docs/phase2/results.json (default: the three candidates)")
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--out", type=Path, default=ROOT / "exports")
    ap.add_argument("--raster-pixels", type=int, default=DEFAULT_RASTER_PIXELS,
                    help="raster size across the pupil diameter (0 = no raster)")
    ap.add_argument("--controller", type=Path, default=None,
                    help="reflex3d_controller checkout: also copy into its config/mask_library/")
    ap.add_argument("--controller-rasters", action="store_true",
                    help="with --controller, copy the .npy rasters too")
    args = ap.parse_args()

    results = json.loads(args.results.read_text())
    system = OpticalSystem()
    library = args.controller / "config" / "mask_library" if args.controller else None
    if library is not None and not (args.controller / "config").is_dir():
        raise SystemExit(f"{args.controller} doesn't look like a reflex3d_controller checkout")

    for config in args.configs:
        if config not in results["configs"]:
            raise SystemExit(f"no configuration {config!r} in {args.results}; "
                             f"available: {', '.join(results['configs'])}")
        res = results["configs"][config]
        family = family_by_name(res["family"], system)
        mask, params = family.mask(res["x"])
        name, label = name_and_label(config)
        path = export_mask(
            args.out, name, mask, params, system, label=label,
            description=(f"{family.describe(res['x'])}. Phase 2 optimum: soft-min EE400 "
                         f"{res['ee_softmin']:.3f}, EE floor {res['ee_min']:.3f}, "
                         f"min Strehl {res['strehl_min']:.3f} over |z| <= 3 um."),
            raster_pixels=args.raster_pixels or None,
            source={"repository": "PhaseMaskOptimization", "results": "docs/phase2/results.json",
                    "configuration": config, "x": res["x"],
                    "metrics": {k: res[k] for k in METRICS if k in res}},
        )
        print(f"{config}: {path.relative_to(ROOT) if path.is_relative_to(ROOT) else path}"
              + (f" + {name}.npy" if args.raster_pixels else ""))
        if library is not None:
            library.mkdir(exist_ok=True)
            shutil.copyfile(path, library / path.name)
            if args.controller_rasters and args.raster_pixels:
                shutil.copyfile(path.with_suffix(".npy"), library / f"{name}.npy")
            print(f"  -> {library / path.name}")


if __name__ == "__main__":
    main()
