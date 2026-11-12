stage 1: 141 CMA-ES runs, 4 workers
  [stage 1] 10/141 runs, 119097 evals, 371 s
  [stage 1] 20/141 runs, 239313 evals, 654 s
  [stage 1] 30/141 runs, 314529 evals, 851 s
  [stage 1] 40/141 runs, 406932 evals, 1110 s
  [stage 1] 50/141 runs, 516549 evals, 1409 s
  [stage 1] 60/141 runs, 599370 evals, 1588 s
  [stage 1] 70/141 runs, 670260 evals, 1755 s
  [stage 1] 80/141 runs, 720320 evals, 1892 s
  [stage 1] 90/141 runs, 779718 evals, 2064 s
  [stage 1] 100/141 runs, 839848 evals, 2221 s
  [stage 1] 110/141 runs, 881933 evals, 2322 s
  [stage 1] 120/141 runs, 922013 evals, 2441 s
  [stage 1] 130/141 runs, 957602 evals, 2525 s
  [stage 1] 140/141 runs, 992692 evals, 2621 s
  [stage 1] 141/141 runs, 996201 evals, 2634 s
  [polish] 10/12 runs, 1036001 evals, 2759 s
  [polish] 12/12 runs, 1048003 evals, 2805 s
  [ablation] 8/8 runs, 1069797 evals, 2875 s
search done: 1069797 evaluations, 2875 s

| configuration | dim | score | p_softmin | p_min | ee_softmin | ee_min | window_fraction_min | eta | strehl_min | ee50_max_core_nm | crlb_max_nm | crlb_max_highbg_nm | feasible | found by |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zoned_zernike6_4 | 19 | 0.2192 | 0.2192 | 0.1778 | 0.2270 | 0.1853 | 0.9535 | 0.7980 | 0.0442 | 1538 | 9.9838 | 19.5468 | yes | polish |
| zoned_zernike4_4 | 15 | 0.2183 | 0.2183 | 0.1715 | 0.2274 | 0.1792 | 0.9459 | 0.7978 | 0.0423 | 1636 | 9.4990 | 19.3343 | yes | truncated+polish |
| zoned_focus_4_sweep | 15 | 0.2182 | 0.2182 | 0.1751 | 0.2264 | 0.1822 | 0.9523 | 0.8005 | 0.0461 | 1469 | 8.9140 | 17.8285 | yes | polish |
| zoned_focus_4_sweep (Phase 2) | 15 | 0.2156 | 0.2156 | 0.1664 | 0.2226 | 0.1744 | 0.9543 | 0.8015 | 0.0420 | 1456 | 8.4810 | 16.3166 | yes | – |
| zoned_zernike6_3 | 14 | 0.2150 | 0.2150 | 0.1745 | 0.2165 | 0.1765 | 0.9890 | 0.7933 | 0.0445 | 1698 | 10.6129 | 20.7078 | yes | warm:zoned_focus_3_sweep |
| zoned_focus_4 | 11 | 0.2129 | 0.2129 | 0.1692 | 0.2192 | 0.1753 | 0.9652 | 0.8009 | 0.0436 | 1592 | 8.3412 | 16.7472 | yes | cold-long |
| zoned_zernike2_4 | 11 | 0.2125 | 0.2125 | 0.1680 | 0.2177 | 0.1729 | 0.9704 | 0.8087 | 0.0453 | 1605 | 9.3040 | 16.9477 | yes | truncated+polish |
| zoned_focus_4 (Phase 2) | 11 | 0.2124 | 0.2124 | 0.1751 | 0.2186 | 0.1813 | 0.9661 | 0.7760 | 0.0380 | 1722 | 10.0733 | 20.6373 | yes | – |
| zoned_focus_3_sweep | 11 | 0.2115 | 0.2115 | 0.1689 | 0.2177 | 0.1744 | 0.9630 | 0.7796 | 0.0445 | 1660 | 11.6051 | 21.7974 | yes | cold-long |
| zoned_zernike4_3 | 11 | 0.2108 | 0.2108 | 0.1668 | 0.2149 | 0.1704 | 0.9705 | 0.7837 | 0.0429 | 1622 | 11.9151 | 22.6553 | yes | truncated+polish |
| zoned_focus_3 | 8 | 0.2078 | 0.2078 | 0.1716 | 0.2111 | 0.1747 | 0.9773 | 0.7823 | 0.0391 | 1890 | 10.2043 | 21.8011 | yes | polish |
| zoned_zernike2_3 | 8 | 0.2072 | 0.2072 | 0.1721 | 0.2104 | 0.1754 | 0.9790 | 0.7781 | 0.0391 | 1919 | 10.5030 | 22.5975 | yes | truncated+polish |
| spherical+defocus | 2 | 0.1650 | 0.1650 | 0.1274 | 0.1668 | 0.1301 | 0.9795 | 0.7834 | 0.0301 | 2311 | 24.0693 | 52.4708 | yes | – |
| no mask | 0 | -1.6654 | 0.0346 | 0.0236 | 0.0348 | 0.0238 | 0.9940 | 0.9044 | 0.0045 | 1405 | 82.6820 | 208 | no | – |

  zoned_zernike6_4: zoned_zernike: radii=(0.5165, 0.708, 0.8781), zone_coefficients=({4: -0.0098, 11: -0.01317, 22: -0.005253}, {1: 0.4189, 4: -0.1504, 11: -0.0596, 22: -0.01964}, {1: 0.3489, 4: 0.1955, 11: -0.04054, 22: -0.09408}, {1: 0.1453, 4: 0.2742, 11: 0.1923, 22: -0.04031}), local_coords=1
  zoned_zernike4_4: zoned_zernike: radii=(0.5305, 0.7414, 0.8896), zone_coefficients=({4: -0.004775, 11: -0.02339}, {1: 0.391, 4: -0.2273, 11: -0.0623}, {1: 0.374, 4: 0.1938, 11: -0.1031}, {1: 0.1295, 4: 0.276, 11: 0.1745}), local_coords=1
  zoned_focus_4_sweep: zoned_focus: radii=(0.5243, 0.7385, 0.8814), depths_nm=(-155, -1626, 1268, 1923), pistons_waves=(0, -1, 1.553, 0.2178), sweeps_nm=(1039, 1757, 4000, -4000)
  zoned_focus_4_sweep (Phase 2): zoned_focus: radii=(0.524, 0.7071, 0.8794), depths_nm=(-304.7, -1822, 1444, 1940), pistons_waves=(0, -0.9959, 1.352, 0.1634), sweeps_nm=(2082, 90.58, -695.2, -3989)
  zoned_zernike6_3: zoned_zernike: radii=(0.5518, 0.8192), zone_coefficients=({4: -0.164, 11: -0.1804, 22: -0.06264}, {1: 1.08, 4: 0.3864, 11: -0.01006, 22: 0.09172}, {1: 0.145, 4: -0.546, 11: -0.07523, 22: 0.06994}), local_coords=1
  zoned_focus_4: zoned_focus: radii=(0.5187, 0.7504, 0.8776), depths_nm=(-239.8, -1716, 1110, 1978), pistons_waves=(0, -0.9932, 0.6354, 0.2157)
  zoned_zernike2_4: zoned_zernike: radii=(0.5136, 0.7108, 0.8752), zone_coefficients=({4: -0.04131}, {1: 0.4045, 4: -0.1987}, {1: 0.3472, 4: 0.2026}, {1: 0.1082, 4: 0.3196}), local_coords=1
  zoned_focus_4 (Phase 2): zoned_focus: radii=(0.5272, 0.6859, 0.8631), depths_nm=(-252.8, 1173, 1519, -2297), pistons_waves=(0, 2, 0.2432, -0.9407)
  zoned_focus_3_sweep: zoned_focus: radii=(0.5, 0.7998), depths_nm=(-292.1, 1622, -2174), pistons_waves=(0, 1.193, 0.9968), sweeps_nm=(2960, 1766, 2925)
  zoned_zernike4_3: zoned_zernike: radii=(0.5504, 0.8004), zone_coefficients=({4: -0.1562, 11: -0.1871}, {1: 1.034, 4: 0.3231, 11: -0.05351}, {1: 0.2032, 4: -0.5256, 11: -0.1716}), local_coords=1
  zoned_focus_3: zoned_focus: radii=(0.4832, 0.8439), depths_nm=(-369.1, 1907, -2703), pistons_waves=(0, -0.838, -0.5229)
  zoned_zernike2_3: zoned_zernike: radii=(0.4861, 0.848), zone_coefficients=({4: -0.0418}, {1: 1.079, 4: 0.4922}, {1: 0.05141, 4: -0.5195}), local_coords=1
  spherical+defocus: zernike: coefficients={4: 0.3715, 11: 0.2924}

done in 2910 s; wrote /home/user/PhaseMaskOptimization/docs/zoned_zernike6/emitted
