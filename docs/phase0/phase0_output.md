system: NA 1.1, 525.0 nm, n 1.33, apodization uniform; 2D grid 1024, 59.7 nm/px
analytic budgets: Int S_axial dz = 925.0 nm (uniform-g formula 901.4, paraxial 1154.1); Int P(r_opt) dz = 1677.9 nm
radial sweeps: 8.4 s for 8 masks x 321 planes (+ wide-range budget checks)

## Budget check (nm)

| mask | Int S_ax ±8 | Int S_ax ±60 | Int P(r_opt) ±8 | Int P(r_opt) ±60 | Int S_peak ±8 | Int EE ±8 |
|---|---|---|---|---|---|---|
| none | 914.1 | 923.6 | 1617.0 | 1669.8 | 950.6 | 1625.9 |
| spherical 2 rad RMS | 909.0 | 923.6 | 1586.0 | 1669.7 | 913.4 | 1624.9 |
| zoned +-1.5pi | 913.2 | 923.6 | 1613.9 | 1669.8 | 924.3 | 1625.2 |
| zoned +-0.5pi | 914.0 | 923.6 | 1616.7 | 1669.8 | 929.0 | 1625.8 |
| zoned +-pi | 913.8 | 923.6 | 1615.6 | 1669.8 | 931.1 | 1625.5 |
| annular 0/pi/0 | 882.7 | 919.0 | 1582.1 | 1664.9 | 891.8 | 1650.8 |
| defocus_sweep 2 um | 913.8 | 923.6 | 1617.2 | 1669.8 | 945.2 | 1626.4 |
| random m=0 Zernikes | 743.2 | 923.2 | 918.0 | 1667.0 | 747.1 | 937.2 |
| analytic (infinite z) | 925.0 | 925.0 | 1677.9 | 1677.9 | not invariant | not invariant |

## In-range summary, |z| <= 3 um (radial model)

| mask | S(0) | EE(0) | eta_P (vs budget) | eta_EE (vs ±8) | EE min | EE mean | flatness | mean S_peak |
|---|---|---|---|---|---|---|---|---|
| none | 1.000 | 0.877 | 0.904 | 0.938 | 0.024 | 0.252 | 0.09 | 0.153 |
| spherical 2 rad RMS | 0.101 | 0.124 | 0.616 | 0.641 | 0.031 | 0.173 | 0.18 | 0.123 |
| zoned +-1.5pi | 0.090 | 0.161 | 0.856 | 0.889 | 0.014 | 0.239 | 0.06 | 0.146 |
| zoned +-0.5pi | 0.811 | 0.738 | 0.901 | 0.935 | 0.017 | 0.251 | 0.07 | 0.149 |
| zoned +-pi | 0.405 | 0.429 | 0.888 | 0.922 | 0.014 | 0.248 | 0.06 | 0.149 |
| annular 0/pi/0 | 0.181 | 0.401 | 0.845 | 0.894 | 0.029 | 0.244 | 0.12 | 0.134 |
| defocus_sweep 2 um | 0.284 | 0.610 | 0.903 | 0.936 | 0.027 | 0.252 | 0.11 | 0.152 |
| random m=0 Zernikes | 0.018 | 0.046 | 0.427 | 0.770 | 0.033 | 0.120 | 0.28 | 0.103 |

flat ceiling over ±3 um: P(r_opt) <= 0.280; over ±4 um: 0.210 (fractions of total light)

2D sweeps: 110.0 s for 3 masks x 161 planes

## Radial vs 2D: max |diff| / max(curve)  [all z | |z|<=3 um]

| mask | strehl | strehl_axial | encircled_energy | ee50_nm | power_weighted_m2 | power_ropt |
|---|---|---|---|---|---|---|
| none | 0.08% / 0.08% | 0.08% / 0.08% | 0.16% / 0.16% | 0.19% / 0.51% | 0.17% / 0.14% | 0.15% / 0.15% |
| spherical 2 rad RMS | 0.46% / 0.46% | 0.46% / 0.46% | 0.48% / 0.48% | 0.19% / 0.49% | 1.65% / 1.40% | 0.48% / 0.48% |
| zoned +-1.5pi | 0.35% / 0.35% | 0.35% / 0.35% | 0.45% / 0.45% | 0.30% / 0.84% | 0.77% / 1.45% | 0.44% / 0.44% |

2D budget integrals over the 2D z grid (trapezoid):
  none                     Int S_axial = 914.9 nm, Int P(r_opt) = 1616.1 nm
  spherical 2 rad RMS      Int S_axial = 909.8 nm, Int P(r_opt) = 1585.0 nm
  zoned +-1.5pi            Int S_axial = 914.1 nm, Int P(r_opt) = 1612.9 nm

## Apodization sensitivity (no mask)

| apodization | Int S_ax dz (nm) | Int P(r_opt) dz (nm) | EE400(0) | EE50(0) nm |
|---|---|---|---|---|
| emission | 907 | 1515 | 0.845 | 130 |
| uniform | 925 | 1678 | 0.877 | 128 |
| illumination | 953 | 1849 | 0.905 | 127 |

wrote /home/user/PhaseMaskOptimization/docs/phase0/results.json and figures/
