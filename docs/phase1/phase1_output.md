targets: {'z_range_nm': (-3000.0, 3000.0), 'dz_nm': 50.0, 's_min': 0.03, 'r_max_nm': 1500.0, 'ee50_range_nm': (-2000.0, 2000.0), 'softmin_p': -10.0, 'penalty': 2.0}
P(r_opt) budget 1678 nm -> flat ceiling 0.280
baselines done (55 s)
NSGA-II tradeoff   zernike_m0: 100 front points (215 s, 10721 evals)
NSGA-II tradeoff   zoned2    : 100 front points (372 s, 18721 evals)
NSGA-II thresholds zernike_m0: 100 front points (530 s, 26721 evals)
NSGA-II thresholds zoned2    : 100 front points (685 s, 34721 evals)

## Best attainable softmin EE400 over |z|<=3 um (rows S_min, cols R_max nm)

| S_min \ R_max | 1000 | 1500 | 2000 | 2500 | 3000 | 4000 | ∞ |
|---|---|---|---|---|---|---|---|
| 0.00 | 0.091 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 |
| 0.01 | 0.091 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 |
| 0.02 | – | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 |
| 0.03 | – | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 | 0.185 |
| 0.04 | – | 0.154 | 0.154 | 0.167 | 0.167 | 0.167 | 0.167 |
| 0.05 | – | 0.126 | 0.144 | 0.159 | 0.159 | 0.159 | 0.159 |
| 0.06 | – | 0.126 | 0.144 | 0.152 | 0.152 | 0.152 | 0.152 |

## Baselines and picks (|z|<=3 um; feasibility at provisional targets)

| mask | ee_softmin | ee_min | ee_mean | eta | ee_flatness | strehl_min | strehl_mean | ee50_max_core_nm | crlb_max_nm | feasible |
|---|---|---|---|---|---|---|---|---|---|---|
| none | 0.034 | 0.024 | 0.252 | 0.904 | 0.094 | 0.005 | 0.153 | 1405 | 137 | no |
| spherical 2 rad RMS | 0.049 | 0.031 | 0.173 | 0.616 | 0.180 | 0.044 | 0.123 | 1450 | 21.110 | yes |
| zoned ±0.5π | 0.026 | 0.017 | 0.251 | 0.901 | 0.067 | 0.004 | 0.149 | 1474 | 154 | no |
| zoned ±π | 0.019 | 0.014 | 0.248 | 0.888 | 0.057 | 0.007 | 0.149 | 1482 | 93.990 | no |
| zoned ±1.5π | 0.021 | 0.014 | 0.239 | 0.856 | 0.060 | 0.007 | 0.146 | 1311 | 55.302 | no |
| annular 0/π/0 | 0.042 | 0.029 | 0.244 | 0.845 | 0.119 | 0.009 | 0.134 | 911 | 28.959 | no |
| defocus_sweep 2 um | 0.035 | 0.027 | 0.252 | 0.903 | 0.107 | 0.005 | 0.152 | 1856 | 160 | no |
| spherical 2 rad + defocus (re-centred) | 0.157 | 0.123 | 0.210 | 0.747 | 0.586 | 0.028 | 0.123 | 2438 | 35.235 | no |
| spherical + defocus (both tuned) | 0.164 | 0.132 | 0.217 | 0.773 | 0.607 | 0.031 | 0.129 | 2326 | 36.261 | no |
| max EE_min (trade-off front) | 0.182 | 0.158 | 0.219 | 0.778 | 0.721 | 0.036 | 0.128 | 1685 | 14.466 | no |
| max softmin EE (unconstrained) | 0.183 | 0.153 | 0.215 | 0.766 | 0.713 | 0.034 | 0.127 | 1517 | 14.760 | no |
| best feasible (0.03, 1500 nm), zernike_m0 | 0.185 | 0.149 | 0.217 | 0.774 | 0.687 | 0.036 | 0.129 | 1498 | 14.753 | yes |
  max EE_min (trade-off front): zernike_m0 x = [0.038, -0.328, -0.136, 0.072]
  max softmin EE (unconstrained): zernike_m0 x = [0.032, -0.329, -0.125, 0.087]
  best feasible (0.03, 1500 nm), zernike_m0: zernike_m0 x = [0.032, -0.318, -0.125, 0.085]
  spherical 2 rad + defocus (re-centred): {'coefficients': {4: 0.41999999999999993, 11: 0.3183098861837907}}
  spherical + defocus (both tuned): {'coefficients': {4: np.float64(0.375), 11: np.float64(0.3)}}

Spearman(softmin EE, -worst CRLB) over 114 trade-off front points: 0.98; (EE_min, -worst CRLB): 0.98

35124 evaluations in 708 s; wrote /home/user/PhaseMaskOptimization/docs/phase1
