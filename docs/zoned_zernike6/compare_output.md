  no mask (control)                      | SLM, mask on the pupil + fringing 0.5 px   softmin P 0.0345 (47 s)
  no mask (control)                      | ideal (continuous)                         softmin P 0.0345 (49 s)
  no mask (control)                      | SLM, mask on the pupil (358.7 px)          softmin P 0.0345 (50 s)
  no mask (control)                      | SLM, mask 362 px on a 358.7 px pupil       softmin P 0.0345 (50 s)
  spherical + defocus (optimised)        | ideal (continuous)                         softmin P 0.1648 (97 s)
  spherical + defocus (optimised)        | SLM, mask on the pupil (358.7 px)          softmin P 0.1646 (98 s)
  spherical + defocus (optimised)        | SLM, mask on the pupil + fringing 0.5 px   softmin P 0.1623 (99 s)
  spherical + defocus (optimised)        | SLM, mask 362 px on a 358.7 px pupil       softmin P 0.1563 (100 s)
  zoned defocus N=4 (defocus + piston)   | ideal (continuous)                         softmin P 0.2122 (146 s)
  zoned defocus N=4 (defocus + piston)   | SLM, mask on the pupil (358.7 px)          softmin P 0.2119 (146 s)
  zoned defocus N=4 (defocus + piston)   | SLM, mask 362 px on a 358.7 px pupil       softmin P 0.2086 (148 s)
  zoned defocus N=4 (defocus + piston)   | SLM, mask on the pupil + fringing 0.5 px   softmin P 0.2062 (149 s)
  zoned defocus N=4 + sweep              | ideal (continuous)                         softmin P 0.2153 (196 s)
  zoned defocus N=4 + sweep              | SLM, mask on the pupil (358.7 px)          softmin P 0.2151 (196 s)
  zoned defocus N=4 + sweep              | SLM, mask 362 px on a 358.7 px pupil       softmin P 0.2092 (197 s)
  zoned defocus N=4 + sweep              | SLM, mask on the pupil + fringing 0.5 px   softmin P 0.2135 (200 s)
  zoned Zernike ≤ 6, N=4 (new)           | SLM, mask on the pupil (358.7 px)          softmin P 0.2187 (242 s)
  zoned Zernike ≤ 6, N=4 (new)           | ideal (continuous)                         softmin P 0.2190 (244 s)
  zoned Zernike ≤ 6, N=4 (new)           | SLM, mask 362 px on a 358.7 px pupil       softmin P 0.2137 (246 s)
  zoned Zernike ≤ 6, N=4 (new)           | SLM, mask on the pupil + fringing 0.5 px   softmin P 0.2119 (249 s)

## Metrics over |z| <= 3 um, by SLM configuration

| mask | configuration | soft-min P400 (emitted) | P400 min | soft-min EE400 (window) | EE400 min | light inside 10 µm, min | Strehl min | worst CRLB | mean CRLB | worst CRLB, high bg | EE50 max ±2 µm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| no mask (control) | ideal (continuous) | 0.034 | 0.024 | 0.035 | 0.024 | 0.994 | 0.005 | 81.7 | 28.1 | 206 | 1406 |
| no mask (control) | SLM, mask on the pupil (358.7 px) | 0.034 | 0.024 | 0.035 | 0.024 | 0.994 | 0.005 | 81.8 | 28.1 | 206 | 1406 |
| no mask (control) | SLM, mask 362 px on a 358.7 px pupil | 0.034 | 0.024 | 0.035 | 0.024 | 0.994 | 0.005 | 81.8 | 28.1 | 206 | 1406 |
| no mask (control) | SLM, mask on the pupil + fringing 0.5 px | 0.034 | 0.024 | 0.035 | 0.024 | 0.994 | 0.005 | 81.8 | 28.1 | 206 | 1406 |
| spherical + defocus (optimised) | ideal (continuous) | 0.165 | 0.127 | 0.167 | 0.130 | 0.980 | 0.030 | 24.1 | 7.6 | 52.5 | 2311 |
| spherical + defocus (optimised) | SLM, mask on the pupil (358.7 px) | 0.165 | 0.127 | 0.166 | 0.130 | 0.979 | 0.030 | 24.1 | 7.6 | 52.4 | 2310 |
| spherical + defocus (optimised) | SLM, mask 362 px on a 358.7 px pupil | 0.156 | 0.103 | 0.158 | 0.103 | 0.981 | 0.033 | 24.4 | 7.6 | 52.5 | 2236 |
| spherical + defocus (optimised) | SLM, mask on the pupil + fringing 0.5 px | 0.162 | 0.120 | 0.169 | 0.124 | 0.956 | 0.026 | 23.2 | 7.6 | 51.4 | 2315 |
| zoned defocus N=4 (defocus + piston) | ideal (continuous) | 0.212 | 0.175 | 0.218 | 0.181 | 0.966 | 0.038 | 10.2 | 6.5 | 20.7 | 1722 |
| zoned defocus N=4 (defocus + piston) | SLM, mask on the pupil (358.7 px) | 0.212 | 0.175 | 0.219 | 0.181 | 0.964 | 0.038 | 10.2 | 6.6 | 20.8 | 1724 |
| zoned defocus N=4 (defocus + piston) | SLM, mask 362 px on a 358.7 px pupil | 0.209 | 0.152 | 0.215 | 0.156 | 0.966 | 0.033 | 11.2 | 6.5 | 24.1 | 1700 |
| zoned defocus N=4 (defocus + piston) | SLM, mask on the pupil + fringing 0.5 px | 0.206 | 0.164 | 0.219 | 0.175 | 0.933 | 0.039 | 10.7 | 6.7 | 22.1 | 1739 |
| zoned defocus N=4 + sweep | ideal (continuous) | 0.215 | 0.166 | 0.222 | 0.174 | 0.954 | 0.044 | 8.5 | 6.3 | 16.5 | 1455 |
| zoned defocus N=4 + sweep | SLM, mask on the pupil (358.7 px) | 0.215 | 0.167 | 0.222 | 0.175 | 0.953 | 0.044 | 8.5 | 6.3 | 16.6 | 1456 |
| zoned defocus N=4 + sweep | SLM, mask 362 px on a 358.7 px pupil | 0.209 | 0.146 | 0.216 | 0.150 | 0.958 | 0.043 | 9.4 | 6.2 | 20.2 | 1443 |
| zoned defocus N=4 + sweep | SLM, mask on the pupil + fringing 0.5 px | 0.214 | 0.169 | 0.223 | 0.178 | 0.944 | 0.042 | 8.3 | 6.3 | 17.1 | 1465 |
| zoned Zernike ≤ 6, N=4 (new) | ideal (continuous) | 0.219 | 0.178 | 0.227 | 0.185 | 0.953 | 0.047 | 10.0 | 6.3 | 19.7 | 1535 |
| zoned Zernike ≤ 6, N=4 (new) | SLM, mask on the pupil (358.7 px) | 0.219 | 0.178 | 0.227 | 0.185 | 0.952 | 0.047 | 10.1 | 6.3 | 19.8 | 1544 |
| zoned Zernike ≤ 6, N=4 (new) | SLM, mask 362 px on a 358.7 px pupil | 0.214 | 0.149 | 0.222 | 0.155 | 0.957 | 0.049 | 12.4 | 6.3 | 26.4 | 1556 |
| zoned Zernike ≤ 6, N=4 (new) | SLM, mask on the pupil + fringing 0.5 px | 0.212 | 0.164 | 0.228 | 0.177 | 0.921 | 0.041 | 11.2 | 6.5 | 22.6 | 1733 |

## Key depths, SLM, mask on the pupil (358.7 px) (P400 / peak Strehl / CRLB nm / EE50 nm)

| mask | -3 µm | -2 µm | -1 µm | +0 µm | +1 µm | +2 µm | +3 µm | worst CRLB |
|---|---|---|---|---|---|---|---|---|
| no mask (control) | 0.024 / 0.007 / 67.8 / 2178 | 0.056 / 0.010 / 43.5 / 1406 | 0.210 / 0.036 / 13.0 / 636 | 0.873 / 1.000 / 2.1 / 127 | 0.210 / 0.036 / 13.0 / 636 | 0.056 / 0.010 / 43.5 / 1406 | 0.024 / 0.007 / 67.8 / 2178 | 0.028 / 0.005 / 81.8 / 2105 @ -2.90 |
| spherical + defocus (optimised) | 0.127 / 0.030 / 24.1 / 3053 | 0.223 / 0.072 / 11.5 / 2310 | 0.304 / 0.139 / 6.7 / 1546 | 0.245 / 0.158 / 5.2 / 929 | 0.147 / 0.087 / 4.9 / 816 | 0.217 / 0.249 / 3.7 / 910 | 0.129 / 0.127 / 5.7 / 1186 | 0.127 / 0.030 / 24.1 / 3053 @ -3.00 |
| zoned defocus N=4 (defocus + piston) | 0.176 / 0.056 / 10.2 / 2039 | 0.214 / 0.045 / 7.0 / 1359 | 0.226 / 0.152 / 7.6 / 1465 | 0.203 / 0.111 / 9.1 / 1635 | 0.221 / 0.041 / 9.2 / 1327 | 0.226 / 0.193 / 5.0 / 1107 | 0.175 / 0.121 / 6.0 / 1680 | 0.176 / 0.056 / 10.2 / 2039 @ -3.00 |
| zoned defocus N=4 + sweep | 0.167 / 0.108 / 8.5 / 2131 | 0.253 / 0.130 / 6.6 / 1456 | 0.204 / 0.109 / 7.2 / 1157 | 0.226 / 0.120 / 8.1 / 1440 | 0.209 / 0.133 / 6.6 / 1326 | 0.227 / 0.102 / 5.2 / 980 | 0.175 / 0.066 / 8.1 / 1370 | 0.167 / 0.108 / 8.5 / 2131 @ -3.00 |
| zoned Zernike ≤ 6, N=4 (new) | 0.178 / 0.103 / 8.8 / 2080 | 0.233 / 0.104 / 7.1 / 1435 | 0.209 / 0.145 / 6.4 / 1145 | 0.229 / 0.108 / 7.3 / 1370 | 0.221 / 0.128 / 6.4 / 1323 | 0.225 / 0.057 / 5.7 / 1544 | 0.178 / 0.078 / 10.1 / 1860 | 0.178 / 0.078 / 10.1 / 1860 @ +3.00 |

## SLM indicators (mask radius 358.7 px)

| mask | phase jumps after the 2π cut (waves) | steepest slope (waves / px) | narrowest zone (px) |
|---|---|---|---|
| spherical + defocus (optimised) | – | 0.029 | – |
| zoned defocus N=4 (defocus + piston) | 0.36, 0.40, 0.08 | 0.020 | 49 |
| zoned defocus N=4 + sweep | 0.31, 0.06, 0.33 | 0.034 | 43 |
| zoned Zernike ≤ 6, N=4 (new) | 0.44, 0.02, 0.31 | 0.030 | 44 |

done in 274 s; wrote /home/user/PhaseMaskOptimization/docs/zoned_zernike6
