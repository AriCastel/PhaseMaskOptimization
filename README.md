# PhaseMaskOptimization

PSF simulation for designing a radially symmetric, extended-depth-of-field
(EDOF) SLM phase mask for single-particle tracking (NA 1.1, 525 nm, water,
target z in [-3, +3] um).

## Layout

```
psfsim/
  config.py      OpticalSystem (NA, wavelength, n, grid, camera, SLM)        [reconstructed]
  zernike.py     Zernike polynomials, Noll indexing                          [reconstructed]
  pupil.py       pupil grid, apodization, exact / paraxial defocus           [reconstructed]
  slm.py         SLM pixelation, wrapping, quantisation                      [reconstructed]
  detector.py    camera area-integration, cropping                           [reconstructed]
  masks.py       registered phase masks (was FourierPlanemasks.py), incl. zoned_focus
  simulator.py   2D FFT forward model, PSFSimulator (was PSFsimulator.py)
  metrics.py     Strehl / EE(r_opt) / PW-M2 per plane, + range summaries (eta, flatness)
  radial.py      Phase 0: fast 1D Debye model for radially symmetric masks + depth budgets
  objective.py   Phase 1: soft-min EE objective, Strehl / EE50 constraints, Evaluator
  crlb.py        Phase 1: lateral CRLB of the camera-sampled PSF (ReflEx3D photometry)
  families.py    Phase 2: parameterised mask families for the optimisers
scripts/phase0_check.py   radial vs 2D cross-check and budget check
scripts/phase1_pareto.py  NSGA-II Pareto study and threshold map
scripts/phase2_zoned.py   zoned_focus search (grid, CMA-ES, Optuna) vs re-optimised baselines
tests/                    Phase 0/1 checks
docs/phase0_report.md     Phase 0 results
docs/phase1_report.md     Phase 1 results (objective, thresholds)
docs/phase2_report.md     Phase 2 results (zoned family vs baselines)
```

Modules marked *reconstructed* were imported by the uploaded files but missing
from the upload. They were rebuilt from how the other modules use them. If you
have the lab versions, drop them in and re-run `scripts/phase0_check.py` and
`pytest`.

## Quick start

```bash
pip install numpy scipy matplotlib pytest pymoo cma optuna
```

```python
import numpy as np
from psfsim import OpticalSystem, RadialPSF

rad = RadialPSF(OpticalSystem(), r_opt_nm=400)
z = np.arange(-8000, 8001, 50.0)
st = rad.simulate(z, "spherical", {"strength_waves": 1 / np.pi})   # 2 rad RMS
st["encircled_energy"], st["strehl"], st["ee50_nm"], st["power_weighted_m2"]
rad.budgets()            # mask-independent Integral S_axial dz and Integral P(r_opt) dz
```

`RadialPSF.simulate` takes the same `mask` / `mask_params` as
`PSFSimulator.simulate` (radially symmetric masks only), plus `phase_waves`
(an array on `rad.rho` or a callable of rho) for free-form phases. A 321-plane sweep
with all metrics takes ~0.3 s (vs ~75 s in 2D); on-axis Strehl alone ~20 ms.

```bash
python scripts/phase0_check.py      # ~2 min; writes docs/phase0/
```

```bash
python -m pytest -q
```
