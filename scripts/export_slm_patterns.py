"""
Export the candidate masks as 8-bit SLM patterns drawn at a given radius.

    python scripts/export_slm_patterns.py --radius 358.7 [--frame 724]

rho = 1 of each mask is placed at `--radius` SLM pixels (use the measured pupil
radius on the SLM, see docs/phase4a_report.md section 2). The pattern is cut
mod 2 pi and quantised to 8 bits (gray = phase * 256 / 2 pi). It is centred,
with the optical axis on the corner of the four central pixels, in a square
frame of `--frame` pixels (default: just large enough). Writes PNG and NPY to
docs/phase4a/slm_patterns/r<radius>/.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from PIL import Image  # noqa: E402

from psfsim import OpticalSystem  # noqa: E402
from psfsim.slm_psf import SLMDevice, SLMPSF  # noqa: E402
from phase4a_slm_compare import candidate_masks  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, required=True, help="mask radius (rho = 1) in SLM pixels")
    ap.add_argument("--frame", type=int, default=None, help="output frame size in pixels (even)")
    a = ap.parse_args()
    sim = SLMPSF(OpticalSystem(), SLMDevice(mask_radius_px=a.radius, pupil_radius_px=a.radius), subsample=1)
    out = ROOT / "docs" / "phase4a" / "slm_patterns" / f"r{a.radius:g}"
    out.mkdir(parents=True, exist_ok=True)
    for name, (m, p) in candidate_masks().items():
        lv = sim.slm_levels(m, p).astype(np.uint8)
        if a.frame:
            if a.frame % 2 or a.frame < lv.shape[0]:
                raise SystemExit("--frame must be even and at least the pattern size")
            pad = (a.frame - lv.shape[0]) // 2
            lv = np.pad(lv, pad)
        fn = name.split(" (")[0].replace(" ", "_").replace("=", "").replace("+", "plus")
        Image.fromarray(lv, mode="L").save(out / f"{fn}.png")  # 8-bit grayscale
        np.save(out / f"{fn}.npy", lv)
        print(f"{out / fn}.png  {lv.shape}  levels {lv.min()}-{lv.max()}")


if __name__ == "__main__":
    main()
