SLM: mask radius 362 px, pupil 358.7 px, 8 bit, 2pi cut; camera 98.5 nm; z -4000..4000 nm step 50
  no mask (control)                  | ideal (continuous)                     softmin 0.0347 EEmin 0.0238 Smin 0.0046 CRLBmax 81.75 nm (46 s)
  no mask (control)                  | SLM, pupil = mask (362 px)             softmin 0.0347 EEmin 0.0238 Smin 0.0046 CRLBmax 81.75 nm (89 s)
  no mask (control)                  | SLM, pupil 3.30 mm (358.7 px)          softmin 0.0347 EEmin 0.0238 Smin 0.0046 CRLBmax 81.75 nm (132 s)
  no mask (control)                  | SLM, pupil 3.30 mm + fringing 0.5 px   softmin 0.0347 EEmin 0.0238 Smin 0.0046 CRLBmax 81.75 nm (175 s)
  zoned_focus N=3                    | ideal (continuous)                     softmin 0.2113 EEmin 0.1740 Smin 0.0379 CRLBmax 10.75 nm (218 s)
  zoned_focus N=3                    | SLM, pupil = mask (362 px)             softmin 0.2114 EEmin 0.1744 Smin 0.0378 CRLBmax 10.68 nm (262 s)
  zoned_focus N=3                    | SLM, pupil 3.30 mm (358.7 px)          softmin 0.2080 EEmin 0.1506 Smin 0.0340 CRLBmax 11.55 nm (306 s)
  zoned_focus N=3                    | SLM, pupil 3.30 mm + fringing 0.5 px   softmin 0.2080 EEmin 0.1503 Smin 0.0332 CRLBmax 11.33 nm (350 s)
  zoned_focus N=4 + sweep            | ideal (continuous)                     softmin 0.2223 EEmin 0.1744 Smin 0.0438 CRLBmax 8.52 nm (393 s)
  zoned_focus N=4 + sweep            | SLM, pupil = mask (362 px)             softmin 0.2225 EEmin 0.1748 Smin 0.0439 CRLBmax 8.48 nm (436 s)
  zoned_focus N=4 + sweep            | SLM, pupil 3.30 mm (358.7 px)          softmin 0.2160 EEmin 0.1503 Smin 0.0434 CRLBmax 9.40 nm (479 s)
  zoned_focus N=4 + sweep            | SLM, pupil 3.30 mm + fringing 0.5 px   softmin 0.2170 EEmin 0.1503 Smin 0.0429 CRLBmax 9.50 nm (522 s)
  spherical + defocus (re-centred)   | ideal (continuous)                     softmin 0.1666 EEmin 0.1300 Smin 0.0301 CRLBmax 24.09 nm (565 s)
  spherical + defocus (re-centred)   | SLM, pupil = mask (362 px)             softmin 0.1665 EEmin 0.1301 Smin 0.0305 CRLBmax 24.04 nm (607 s)
  spherical + defocus (re-centred)   | SLM, pupil 3.30 mm (358.7 px)          softmin 0.1577 EEmin 0.1032 Smin 0.0329 CRLBmax 24.41 nm (652 s)
  spherical + defocus (re-centred)   | SLM, pupil 3.30 mm + fringing 0.5 px   softmin 0.1540 EEmin 0.0978 Smin 0.0274 CRLBmax 24.76 nm (696 s)

## Metrics over |z| <= 3 um, by SLM configuration

| mask | configuration | soft-min EE | EE min | EE mean | S min | S mean | worst CRLB | mean CRLB | worst CRLB hi-bg | EE50 max ±2 | S≥0.03 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| no mask (control) | ideal (continuous) | 0.035 | 0.024 | 0.254 | 0.005 | 0.154 | 81.7 | 28.1 | 206 | 1406 | no |
| no mask (control) | SLM, pupil = mask (362 px) | 0.035 | 0.024 | 0.254 | 0.005 | 0.154 | 81.7 | 28.1 | 206 | 1406 | no |
| no mask (control) | SLM, pupil 3.30 mm (358.7 px) | 0.035 | 0.024 | 0.254 | 0.005 | 0.154 | 81.8 | 28.1 | 206 | 1406 | no |
| no mask (control) | SLM, pupil 3.30 mm + fringing 0.5 px | 0.035 | 0.024 | 0.254 | 0.005 | 0.154 | 81.8 | 28.1 | 206 | 1406 | no |
| zoned_focus N=3 | ideal (continuous) | 0.211 | 0.174 | 0.222 | 0.038 | 0.129 | 10.8 | 6.355 | 20.0 | 1682 | yes |
| zoned_focus N=3 | SLM, pupil = mask (362 px) | 0.211 | 0.174 | 0.223 | 0.038 | 0.129 | 10.7 | 6.358 | 19.8 | 1682 | yes |
| zoned_focus N=3 | SLM, pupil 3.30 mm (358.7 px) | 0.208 | 0.151 | 0.225 | 0.034 | 0.130 | 11.6 | 6.323 | 22.0 | 1655 | yes |
| zoned_focus N=3 | SLM, pupil 3.30 mm + fringing 0.5 px | 0.208 | 0.150 | 0.227 | 0.033 | 0.127 | 11.3 | 6.424 | 21.6 | 1688 | yes |
| zoned_focus N=4 + sweep | ideal (continuous) | 0.222 | 0.174 | 0.231 | 0.044 | 0.132 | 8.520 | 6.256 | 16.5 | 1455 | yes |
| zoned_focus N=4 + sweep | SLM, pupil = mask (362 px) | 0.222 | 0.175 | 0.231 | 0.044 | 0.132 | 8.485 | 6.257 | 16.6 | 1456 | yes |
| zoned_focus N=4 + sweep | SLM, pupil 3.30 mm (358.7 px) | 0.216 | 0.150 | 0.234 | 0.043 | 0.134 | 9.397 | 6.181 | 20.2 | 1443 | yes |
| zoned_focus N=4 + sweep | SLM, pupil 3.30 mm + fringing 0.5 px | 0.217 | 0.150 | 0.235 | 0.043 | 0.132 | 9.495 | 6.229 | 20.4 | 1444 | yes |
| spherical + defocus (re-centred) | ideal (continuous) | 0.167 | 0.130 | 0.221 | 0.030 | 0.131 | 24.1 | 7.626 | 52.5 | 2311 | yes |
| spherical + defocus (re-centred) | SLM, pupil = mask (362 px) | 0.166 | 0.130 | 0.221 | 0.030 | 0.132 | 24.0 | 7.621 | 52.4 | 2310 | yes |
| spherical + defocus (re-centred) | SLM, pupil 3.30 mm (358.7 px) | 0.158 | 0.103 | 0.224 | 0.033 | 0.134 | 24.4 | 7.612 | 52.5 | 2236 | yes |
| spherical + defocus (re-centred) | SLM, pupil 3.30 mm + fringing 0.5 px | 0.154 | 0.098 | 0.227 | 0.027 | 0.132 | 24.8 | 7.683 | 54.7 | 2241 | no |

## Change from ideal to real SLM (pupil 3.30 mm), and with fringing

| mask | soft-min EE | EE min | S min | worst CRLB |
|---|---|---|---|---|
| no mask (control) | +0.0% / +0.0% | +0.0% / +0.0% | -0.0% / -0.0% | +0.0% / +0.0% |
| zoned_focus N=3 | -1.6% / -1.6% | -13.5% / -13.6% | -10.3% / -12.6% | +7.5% / +5.4% |
| zoned_focus N=4 + sweep | -2.8% / -2.4% | -13.8% / -13.8% | -0.8% / -2.0% | +10.3% / +11.4% |
| spherical + defocus (re-centred) | -5.3% / -7.6% | -20.6% / -24.7% | +9.4% / -8.8% | +1.3% / +2.8% |

## Key depths on the real SLM (EE400 / peak Strehl / CRLB nm / EE50 nm)

| mask | -4 µm | -3 µm | -2 µm | -1 µm | +0 µm | +1 µm | +2 µm | +3 µm | +4 µm | worst CRLB | best CRLB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| no mask (control) | 0.013 / 0.005 / 98.1 / 2907 | 0.024 / 0.007 / 67.8 / 2178 | 0.056 / 0.010 / 43.5 / 1406 | 0.211 / 0.036 / 13.0 / 636 | 0.877 / 1.000 / 2.1 / 127 | 0.211 / 0.036 / 13.0 / 636 | 0.056 / 0.010 / 43.5 / 1406 | 0.024 / 0.007 / 67.8 / 2178 | 0.013 / 0.005 / 98.1 / 2907 | 0.028 / 0.005 / 81.8 / 2105 @-2.90 | 0.877 / 1.000 / 2.1 / 127 @+0.00 |
| zoned_focus N=3 | 0.021 / 0.011 / 20.7 / 2001 | 0.166 / 0.083 / 11.6 / 1423 | 0.218 / 0.177 / 5.3 / 878 | 0.212 / 0.104 / 5.2 / 868 | 0.232 / 0.138 / 6.6 / 1366 | 0.250 / 0.192 / 5.8 / 1653 | 0.207 / 0.085 / 7.0 / 1625 | 0.151 / 0.034 / 8.1 / 2175 | 0.050 / 0.023 / 22.5 / 2828 | 0.166 / 0.083 / 11.6 / 1423 @-3.00 | 0.234 / 0.083 / 4.6 / 802 @-1.30 |
| zoned_focus N=4 + sweep | 0.016 / 0.009 / 23.6 / 3062 | 0.164 / 0.099 / 8.8 / 2129 | 0.267 / 0.149 / 6.2 / 1443 | 0.215 / 0.124 / 6.6 / 1104 | 0.235 / 0.128 / 7.7 / 1387 | 0.221 / 0.156 / 6.1 / 1309 | 0.224 / 0.132 / 5.2 / 951 | 0.150 / 0.043 / 9.4 / 1467 | 0.019 / 0.006 / 36.3 / 2555 | 0.150 / 0.043 / 9.4 / 1467 @+3.00 | 0.209 / 0.278 / 4.1 / 958 @+2.30 |
| spherical + defocus (re-centred) | 0.069 / 0.014 / 54.4 / 3700 | 0.127 / 0.033 / 24.4 / 2997 | 0.225 / 0.077 / 11.4 / 2236 | 0.311 / 0.152 / 6.4 / 1485 | 0.255 / 0.178 / 5.0 / 884 | 0.148 / 0.084 / 4.7 / 795 | 0.232 / 0.266 / 3.7 / 888 | 0.103 / 0.113 / 6.2 / 1241 | 0.020 / 0.019 / 21.6 / 1874 | 0.127 / 0.033 / 24.4 / 2997 @-3.00 | 0.228 / 0.263 / 3.7 / 881 @+1.95 |

done in 714 s; wrote /home/user/PhaseMaskOptimization/docs/phase4a
