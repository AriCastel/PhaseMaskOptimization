stage 1: 135 CMA-ES runs, 4 workers
  [stage 1] 10/135 runs, 120120 evals, 386 s
  [stage 1] 20/135 runs, 240336 evals, 660 s
  [stage 1] 30/135 runs, 314292 evals, 863 s
  [stage 1] 40/135 runs, 408015 evals, 1104 s
  [stage 1] 50/135 runs, 515025 evals, 1401 s
  [stage 1] 60/135 runs, 598761 evals, 1586 s
  [stage 1] 70/135 runs, 664294 evals, 1784 s
  [stage 1] 80/135 runs, 722943 evals, 1908 s
  [stage 1] 90/135 runs, 783067 evals, 2078 s
  [stage 1] 100/135 runs, 833172 evals, 2191 s
  [stage 1] 110/135 runs, 873252 evals, 2316 s
  [stage 1] 120/135 runs, 910837 evals, 2398 s
  [stage 1] 130/135 runs, 945927 evals, 2507 s
  [stage 1] 135/135 runs, 963472 evals, 2544 s
  [polish] 10/12 runs, 1070362 evals, 2869 s
  [polish] 12/12 runs, 1096992 evals, 2958 s
  [ablation] 8/8 runs, 1135133 evals, 3077 s
search done: 1135133 evaluations, 3077 s

| configuration | dim | score | ee_softmin | ee_min | ee_mean | eta | ee_flatness | strehl_min | ee50_max_core_nm | crlb_max_nm | crlb_max_highbg_nm | feasible | found by |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zoned_zernike4_4 | 15 | 0.2599 | 0.2599 | 0.2098 | 0.2748 | 0.6128 | 0.7636 | 0.0325 | 2574 | 17.2919 | 37.9455 | yes | truncated+polish |
| zoned_zernike6_4 | 19 | 0.2594 | 0.2594 | 0.1989 | 0.2758 | 0.6529 | 0.7212 | 0.0361 | 2649 | 18.7104 | 37.9474 | yes | polish |
| zoned_zernike2_4 | 11 | 0.2522 | 0.2522 | 0.2091 | 0.2644 | 0.6002 | 0.7907 | 0.0301 | 3085 | 17.2939 | 35.7325 | yes | truncated+polish |
| zoned_zernike6_3 | 14 | 0.2427 | 0.2427 | 0.1961 | 0.2517 | 0.7022 | 0.7789 | 0.0383 | 2712 | 13.9234 | 27.0438 | yes | warm:zoned_focus_3_sweep |
| zoned_zernike4_3 | 11 | 0.2308 | 0.2308 | 0.1817 | 0.2414 | 0.7178 | 0.7530 | 0.0313 | 2740 | 13.3980 | 28.9419 | yes | truncated+polish |
| zoned_focus_4_sweep | 15 | 0.2265 | 0.2265 | 0.1822 | 0.2320 | 0.8000 | 0.7853 | 0.0469 | 1477 | 9.0313 | 17.9471 | yes | polish |
| zoned_zernike2_3 | 8 | 0.2237 | 0.2237 | 0.1701 | 0.2587 | 0.6923 | 0.6576 | 0.0326 | 2908 | 18.7814 | 39.6366 | yes | truncated+polish |
| zoned_focus_4_sweep (Phase 2) | 15 | 0.2226 | 0.2226 | 0.1744 | 0.2309 | 0.8015 | 0.7552 | 0.0420 | 1456 | 8.4810 | 16.3166 | yes | – |
| zoned_focus_4 | 11 | 0.2186 | 0.2186 | 0.1805 | 0.2232 | 0.7768 | 0.8089 | 0.0380 | 1721 | 10.0012 | 20.5703 | yes | polish |
| zoned_focus_4 (Phase 2) | 11 | 0.2186 | 0.2186 | 0.1813 | 0.2232 | 0.7760 | 0.8121 | 0.0380 | 1722 | 10.0733 | 20.6373 | yes | – |
| zoned_focus_3_sweep | 11 | 0.2180 | 0.2180 | 0.1739 | 0.2242 | 0.7778 | 0.7759 | 0.0446 | 1663 | 11.5869 | 21.8193 | yes | cold-long |
| zoned_focus_3 | 8 | 0.2115 | 0.2115 | 0.1740 | 0.2224 | 0.7751 | 0.7826 | 0.0378 | 1683 | 10.7321 | 19.9308 | yes | phase2 |
| spherical+defocus | 2 | 0.1668 | 0.1668 | 0.1301 | 0.2206 | 0.7834 | 0.5895 | 0.0301 | 2311 | 24.0693 | 52.4708 | yes | – |
| no mask | 0 | -1.6651 | 0.0348 | 0.0238 | 0.2537 | 0.9044 | 0.0937 | 0.0045 | 1405 | 82.6820 | 208 | no | – |

  zoned_zernike4_4: zoned_zernike: radii=(0.4928, 0.7922, 0.9135), zone_coefficients=({4: 0.2129, 11: 0.1463}, {1: 0.003931, 4: -0.3153, 11: 0.04827}, {1: 1.225, 4: 1.942, 11: -0.08084}, {1: 1.701, 4: 1.627, 11: 0.1303}), local_coords=1
  zoned_zernike6_4: zoned_zernike: radii=(0.5053, 0.7772, 0.9017), zone_coefficients=({4: 0.04526, 11: 0.01841, 22: 0.0009864}, {1: 0.2729, 4: -0.2575, 11: 0.006341, 22: 0.1634}, {1: 2, 4: 2, 11: -0.07005, 22: 0.1428}, {1: 1.57, 4: 1.077, 11: 0.4798, 22: -0.3225}), local_coords=1
  zoned_zernike2_4: zoned_zernike: radii=(0.5263, 0.8001, 0.931), zone_coefficients=({4: 0.1941}, {1: -0.4164, 4: -0.4248}, {1: 1.618, 4: 2}, {1: 1.799, 4: 1.171}), local_coords=1
  zoned_zernike6_3: zoned_zernike: radii=(0.6, 0.8635), zone_coefficients=({4: -0.224, 11: -0.04046, 22: 0.09406}, {1: 1.404, 4: 0.3515, 11: 0.04661, 22: -0.08124}, {1: 0.9148, 4: -1.963, 11: 1, 22: 0.5}), local_coords=1
  zoned_zernike4_3: zoned_zernike: radii=(0.5842, 0.8697), zone_coefficients=({4: -0.2619, 11: -0.096}, {1: 1.375, 4: 0.4169, 11: 0.04772}, {1: 0.8193, 4: -1.97, 11: 1}), local_coords=1
  zoned_focus_4_sweep: zoned_focus: radii=(0.5234, 0.7392, 0.8814), depths_nm=(-137.6, -1603, 1200, 1934), pistons_waves=(0, -1, 1.605, 0.2149), sweeps_nm=(721.7, 2192, 4000, -4000)
  zoned_zernike2_3: zoned_zernike: radii=(0.5616, 0.8686), zone_coefficients=({4: -0.2137}, {1: 1.472, 4: 0.4774}, {1: 0.609, 4: -2}), local_coords=1
  zoned_focus_4_sweep (Phase 2): zoned_focus: radii=(0.524, 0.7071, 0.8794), depths_nm=(-304.7, -1822, 1444, 1940), pistons_waves=(0, -0.9959, 1.352, 0.1634), sweeps_nm=(2082, 90.58, -695.2, -3989)
  zoned_focus_4: zoned_focus: radii=(0.5279, 0.6864, 0.8531), depths_nm=(-248.1, 1177, 1561, -2314), pistons_waves=(0, 2, 1.221, -0.9222)
  zoned_focus_4 (Phase 2): zoned_focus: radii=(0.5272, 0.6859, 0.8631), depths_nm=(-252.8, 1173, 1519, -2297), pistons_waves=(0, 2, 0.2432, -0.9407)
  zoned_focus_3_sweep: zoned_focus: radii=(0.5048, 0.7991), depths_nm=(-224.8, 1587, -2161), pistons_waves=(0, -0.7894, 0.002402), sweeps_nm=(2405, 2273, 3064)
  zoned_focus_3: zoned_focus: radii=(0.4812, 0.813), depths_nm=(164.7, -1592, 2351), pistons_waves=(0, 1.748, 0.8411)
  spherical+defocus: zernike: coefficients={4: 0.3715, 11: 0.2924}

done in 3111 s; wrote /home/user/PhaseMaskOptimization/docs/zoned_zernike6
