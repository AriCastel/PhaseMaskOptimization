305 tasks (114376 grid points), 4 workers
  20/305 tasks, 16568 evals, 75 s
  40/305 tasks, 34976 evals, 132 s
  60/305 tasks, 42956 evals, 152 s
  80/305 tasks, 48068 evals, 169 s
  100/305 tasks, 61738 evals, 210 s
  120/305 tasks, 80786 evals, 268 s
  140/305 tasks, 107971 evals, 344 s
  160/305 tasks, 144672 evals, 448 s
  180/305 tasks, 189639 evals, 596 s
  200/305 tasks, 244859 evals, 773 s
  220/305 tasks, 307044 evals, 964 s
  240/305 tasks, 380435 evals, 1172 s
  260/305 tasks, 454108 evals, 1380 s
  280/305 tasks, 531237 evals, 1588 s
  300/305 tasks, 602718 evals, 1777 s
  305/305 tasks, 620588 evals, 1828 s
search done: 645335 evaluations, 1900 s

## Best mask per configuration (reference evaluator, |z| <= 3 um)

| configuration | dim | score | ee_softmin | ee_min | ee_mean | eta | ee_flatness | strehl_min | strehl_mean | ee50_max_core_nm | crlb_max_nm | crlb_max_highbg_nm | feasible | CMA hit rate | Optuna |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zoned_focus_4_sweep | 15 | 0.223 | 0.223 | 0.174 | 0.231 | 0.802 | 0.755 | 0.042 | 0.132 | 1456 | 8.481 | 16.317 | yes | 0.05 | – |
| zoned_focus_5_sweep | 19 | 0.221 | 0.221 | 0.174 | 0.233 | 0.763 | 0.748 | 0.031 | 0.123 | 2291 | 11.481 | 24.755 | yes | 0.10 | – |
| zoned_focus_5 | 14 | 0.221 | 0.221 | 0.177 | 0.225 | 0.768 | 0.786 | 0.038 | 0.123 | 1850 | 10.245 | 21.667 | yes | 0.05 | – |
| zoned_focus_4 | 11 | 0.219 | 0.219 | 0.181 | 0.223 | 0.776 | 0.812 | 0.038 | 0.128 | 1722 | 10.073 | 20.637 | yes | 0.10 | – |
| zoned_focus_3_sweep | 11 | 0.218 | 0.218 | 0.174 | 0.224 | 0.777 | 0.776 | 0.046 | 0.128 | 1655 | 11.703 | 22.076 | yes | 0.10 | 0.210 |
| annular_binary_3 | 4 | 0.214 | 0.214 | 0.175 | 0.224 | 0.765 | 0.781 | 0.042 | 0.128 | 1814 | 8.261 | 17.470 | yes | 0.80 | – |
| zoned_focus_3 | 8 | 0.212 | 0.212 | 0.174 | 0.222 | 0.775 | 0.783 | 0.038 | 0.129 | 1683 | 10.732 | 19.931 | yes | 0.10 | 0.196 |
| zoned_focus_2_sweep | 7 | 0.207 | 0.207 | 0.163 | 0.218 | 0.764 | 0.749 | 0.041 | 0.120 | 2711 | 12.679 | 25.145 | yes | 0.35 | – |
| zernike_m0 | 4 | 0.185 | 0.185 | 0.147 | 0.201 | 0.711 | 0.730 | 0.030 | 0.122 | 2484 | 19.510 | 42.566 | yes | 0.40 | – |
| annular_binary_2 | 3 | 0.181 | 0.181 | 0.135 | 0.218 | 0.758 | 0.620 | 0.060 | 0.128 | 2008 | 13.688 | 27.161 | yes | 0.20 | – |
| layer_cake_5 | 4 | 0.180 | 0.180 | 0.130 | 0.233 | 0.799 | 0.560 | 0.043 | 0.132 | 2500 | 17.661 | 36.736 | yes | 1.00 | – |
| spherical+defocus | 2 | 0.167 | 0.167 | 0.130 | 0.221 | 0.783 | 0.590 | 0.030 | 0.131 | 2311 | 24.069 | 52.471 | yes | 1.00 | – |
| zoned_focus_2 | 5 | 0.166 | 0.166 | 0.127 | 0.228 | 0.809 | 0.556 | 0.033 | 0.137 | 3232 | 19.973 | 42.728 | yes | 0.55 | 0.164 |
| zoned_zernike2 | 4 | 0.131 | 0.164 | 0.125 | 0.227 | 0.804 | 0.553 | 0.030 | 0.136 | 3316 | 20.781 | 44.851 | no | 0.80 | – |
| annular_binary_1 | 2 | 0.004 | 0.126 | 0.095 | 0.232 | 0.817 | 0.410 | 0.028 | 0.137 | 1618 | 22.992 | 52.196 | no | 0.90 | – |
| layer_cake_4 | 3 | -0.706 | 0.056 | 0.046 | 0.247 | 0.863 | 0.185 | 0.019 | 0.145 | 1088 | 33.030 | 77.784 | no | 0.00 | – |
| layer_cake_2 | 1 | -0.994 | 0.056 | 0.047 | 0.246 | 0.866 | 0.189 | 0.014 | 0.150 | 1056 | 29.064 | 68.027 | no | 0.00 | – |
| layer_cake_3 | 2 | -1.227 | 0.022 | 0.017 | 0.239 | 0.826 | 0.072 | 0.011 | 0.139 | 1352 | 16.199 | 38.505 | no | 0.00 | – |

  zoned_focus_4_sweep: zoned_focus: radii=(0.524, 0.7071, 0.8794), depths_nm=(-304.7, -1822, 1444, 1940), pistons_waves=(0, -0.9959, 1.352, 0.1634), sweeps_nm=(2082, 90.58, -695.2, -3989)
  zoned_focus_5_sweep: zoned_focus: radii=(0.5859, 0.8669, 0.8805, 0.8955), depths_nm=(1455, -1882, 3791, -1381, 2510), pistons_waves=(0, 1.923, 0.9307, 0.4421, -0.05638), sweeps_nm=(-3891, 641.5, 680.3, -2409, 3908)
  zoned_focus_5: zoned_focus: radii=(0.5205, 0.7769, 0.8766, 0.9699), depths_nm=(200.1, -1504, 3422, -25.72, 2471), pistons_waves=(0, 0.9392, 2, 0.9884, 0.7774)
  zoned_focus_4: zoned_focus: radii=(0.5272, 0.6859, 0.8631), depths_nm=(-252.8, 1173, 1519, -2297), pistons_waves=(0, 2, 0.2432, -0.9407)
  zoned_focus_3_sweep: zoned_focus: radii=(0.5031, 0.7968), depths_nm=(-232.4, 1581, -2129), pistons_waves=(0, 0.2144, 1.965), sweeps_nm=(2495, 2387, 3190)
  annular_binary_3: annular: radii=(0.5206, 0.6591, 0.7861), phases_waves=(0, 0.5, 0, 0.5)
  zoned_focus_3: zoned_focus: radii=(0.4812, 0.813), depths_nm=(164.7, -1592, 2351), pistons_waves=(0, 1.748, 0.8411)
  zoned_focus_2_sweep: zoned_focus: radii=(0.5875), depths_nm=(-1565, 2290), pistons_waves=(0, -0.2886), sweeps_nm=(1697, -4000)
  zernike_m0: zernike: coefficients={4: -0.4358, 11: -0.3076, 22: 0.002439, 37: -0.02533}
  annular_binary_2: annular: radii=(0.4396, 0.6485), phases_waves=(0, 0.5, 0)
  layer_cake_5: zoned_focus: radii=(0.4948, 0.6836, 0.8169, 0.9192), depths_nm=(0, 0, 0, 0, 0), pistons_waves=(0, 1.587, 0.1075, 0.4233, -0.4987)
  spherical+defocus: zernike: coefficients={4: 0.3715, 11: 0.2924}
  zoned_focus_2: zoned_focus: radii=(0.5485), depths_nm=(1307, -2160), pistons_waves=(0, -0.2783)
  zoned_zernike2: zoned_zernike: radii=(0.5507), zone_coefficients=({4: 0.5896}, {4: -1.258, 1: 0.06103})
  annular_binary_1: annular: radii=(0.4622), phases_waves=(0, 0.5)
  layer_cake_4: zoned_focus: radii=(0.5501, 0.755, 0.8959), depths_nm=(0, 0, 0, 0), pistons_waves=(0, 1.972, 0.4703, 0.5982)
  layer_cake_2: zoned_focus: radii=(0.755), depths_nm=(0, 0), pistons_waves=(0, 0.4985)
  layer_cake_3: zoned_focus: radii=(0.629, 0.8538), depths_nm=(0, 0, 0), pistons_waves=(0, 1.5, 1)

done in 1940 s; wrote /home/user/PhaseMaskOptimization/docs/phase2
exit 0
