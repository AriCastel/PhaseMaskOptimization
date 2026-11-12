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
  export.py      mask files for the microscope: parameters (.json) + rasterised matrix (.npy)
scripts/phase0_check.py   radial vs 2D cross-check and budget check
scripts/phase1_pareto.py  NSGA-II Pareto study and threshold map
scripts/phase2_zoned.py   zoned_focus search (grid, CMA-ES, Optuna) vs re-optimised baselines
scripts/export_masks.py   export optimised masks for the microscope (reflex3d_controller)
exports/                  the exported candidate masks (.json; the .npy rasters are generated)
tests/                    Phase 0/1 checks
docs/phase0_report.md     Phase 0 results
docs/phase1_report.md     Phase 1 results (objective, thresholds)
docs/phase2_report.md     Phase 2 results (zoned family vs baselines)
docs/phase4a_report.md    Phase 4a: candidates on the real SLM
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

## Exporting masks to the microscope

```bash
python scripts/export_masks.py
```

writes the three masks carried to the bench into `exports/`, each as
`<name>.json` (parameters and provenance) plus `<name>.npy` (rasterised matrix):

| name | mask | from `docs/phase2/results.json` |
|---|---|---|
| `zoned_focus_3` | 3-zone `zoned_focus` | `zoned_focus_3` |
| `zoned_focus_4_sweep` | 4-zone `zoned_focus` with per-zone sweep | `zoned_focus_4_sweep` |
| `spherical_defocus` | re-centred spherical, Z11 + Z4 | `spherical+defocus` |

`--configs` picks other configurations (e.g. `annular_binary_3`), and

```bash
python scripts/export_masks.py --controller ../reflex3d_controller
```

also copies the `.json` files into that checkout's `config/mask_library/`, where
the controller lists them in its 'SLM phase mask' panel (add
`--controller-rasters` to copy the matrices too).

The two forms, both in waves and both described in `psfsim/export.py`:

- **Parameters** (`"parametric"` in the `.json`): `zoned_focus` (radii, depths,
  pistons, sweeps, with the NA / n / wavelength they were optimised for),
  `zernike` (explicit `(n, m)` modes in waves RMS) or `annular` (radii, phases).
  The controller evaluates these analytically on each channel's pupil, so they
  are exact for any pupil size and centre on the SLM. This is what it uses when
  a file has both forms.
- **Rasterised matrix** (`<name>.npy`, `np.load` gives it): the unwrapped phase
  as a square float32 array across the pupil diameter (1024 px by default,
  `--raster-pixels`), 0 outside the pupil; element `[i, j]` is at
  `x = (2j + 1)/N - 1`, `y = (2i + 1)/N - 1`. Any registered mask can be
  exported this way, including ones with no parametric form in the controller
  (`cubic`, `zoned_zernike`, future free-form phases), which are then
  raster-only. The `.npy` files are not committed (4 MB each); the script
  regenerates them.

From code, for any mask:

```python
from psfsim.export import export_mask
export_mask("exports", "my_mask", "zoned_focus", {"radii": (0.5, 0.8), "depths_nm": (0, -1500, 2000)})
```

