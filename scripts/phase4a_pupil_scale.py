"""
Phase 4a: sensitivity to the pupil / drawn-mask radius ratio k = R_pupil / R_mask.

The mask is designed in normalised pupil radius rho. Drawn on the SLM with
rho = 1 at R_mask pixels while the true pupil edge is at R_pupil, the light at
pupil radius rho_p sees the mask at rho_mask = k * rho_p (k = 358.7 / 362 =
0.991 for the 3.30 mm pupil and the 362 px mask). Radial model, reference
resolution. Writes docs/phase4a/pupil_scale.json and figures/pupil_scale.png.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from psfsim import OpticalSystem, RadialPSF  # noqa: E402
from psfsim.masks import build_mask  # noqa: E402
from psfsim.objective import Evaluator, Targets  # noqa: E402
from psfsim.pupil import radial_pupil_grid  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from phase4a_slm_compare import candidate_masks  # noqa: E402

OUT = ROOT / "docs" / "phase4a"
K = np.round(np.arange(0.95, 1.0501, 0.005), 4)
K_REAL = 3300.0 / 9.2 / 362.0


def main():
    ev = Evaluator(RadialPSF(OpticalSystem()), Targets(dz_nm=25.0))
    res = {}
    for name, (m, p) in candidate_masks().items():
        if m is None:
            continue
        rows = []
        for k in list(K) + [K_REAL]:
            ph = build_mask(m, radial_pupil_grid(ev.rad.rho * k), ev.rad.system, p)
            r = ev.evaluate(ph, with_crlb=True)
            rows.append({"k": float(k), **{c: float(r[c]) for c in
                         ("ee_softmin", "ee_min", "strehl_min", "crlb_max_nm", "eta")}})
        res[name] = rows
        real = rows[-1]
        print(f"{name:34s} k=1: softmin {rows[10]['ee_softmin']:.3f} EEmin {rows[10]['ee_min']:.3f} "
              f"CRLBmax {rows[10]['crlb_max_nm']:.2f} | k={K_REAL:.4f}: softmin {real['ee_softmin']:.3f} "
              f"EEmin {real['ee_min']:.3f} CRLBmax {real['crlb_max_nm']:.2f}")
    (OUT / "pupil_scale.json").write_text(json.dumps({"k_real": K_REAL, "results": res}, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"zoned_focus N=3": "#2f855a", "zoned_focus N=4 + sweep": "#2b6cb0",
              "spherical + defocus (re-centred)": "#c05621"}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for name, rows in res.items():
        rows = rows[:-1]
        k = [r["k"] for r in rows]
        for ax, key in zip(axes, ("ee_softmin", "ee_min", "crlb_max_nm")):
            ax.plot(k, [r[key] for r in rows], "o-", ms=3, color=colors[name], label=name)
    for ax, lab in zip(axes, ("soft-min EE400 (±3 µm)", "EE400 min (±3 µm)", "worst lateral CRLB (nm, ±3 µm)")):
        ax.axvline(1.0, color="k", lw=0.8, ls=":")
        ax.axvline(K_REAL, color="#e53e3e", lw=1.0, ls="--")
        ax.set_xlabel("k = pupil radius / drawn mask radius")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[0].text(K_REAL, axes[0].get_ylim()[1], "3.30 mm / 362 px ", color="#e53e3e", fontsize=8, rotation=90, ha="right",
                 va="top")
    axes[2].legend(fontsize=8, loc="upper right")
    fig.suptitle("Sensitivity to the pupil / mask size ratio (k < 1: mask drawn larger than the pupil)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "pupil_scale.png", dpi=140)


if __name__ == "__main__":
    main()
