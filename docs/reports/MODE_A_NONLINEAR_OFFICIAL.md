# Mode A step 2: non-linear probes on TEDDY's pooled layers

Written by `scripts/mode_a_nonlinear_probes.py` (mode_a_nonlinear_probes v1), 2026-09-30 17:28:46. Status: **complete**.

Same cells, targets, splits, look-alike pairs (all 1506, no_gdT158 453), bootstrap and head as `README.md` (the linear run); checks: pairs identical to the linear run = True, head gap ratios reproduced (max abs diff 0.0000).

MLP: 2 x 256 ReLU, dropout 0.1, AdamW lr 0.001, early stopping on val (patience 8, max 100 epochs), seeds [0, 1, 2]; 'MLP' rows = mean of the seeds' predictions, per-seed ratios in the JSON. kNN: k per protein on val from [5, 10, 20, 50, 100, 200]; cos_std = cosine on train-standardised features, cos_raw = cosine on raw features (at z: the pair-selection geometry).

## site4 Pearson (p95-scaled, clipped)

| depth | read-out | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|---|
| z | linear (ridge) | 0.718 | 0.687 | 0.582 | 0.896 |
| z | MLP | 0.779 | 0.712 | 0.613 | 0.932 |
| z | kNN cos_std | 0.706 | 0.674 | 0.592 | 0.902 |
| z | kNN cos_raw | 0.708 | 0.676 | 0.591 | 0.903 |
| 12 | linear (ridge) | 0.717 | 0.685 | 0.582 | 0.896 |
| 12 | MLP | 0.782 | 0.710 | 0.615 | 0.928 |
| 12 | kNN cos_std | 0.706 | 0.675 | 0.592 | 0.902 |
| 12 | kNN cos_raw | 0.708 | 0.676 | 0.591 | 0.903 |
| 9 | linear (ridge) | 0.659 | 0.657 | 0.507 | 0.846 |
| 9 | MLP | 0.724 | 0.688 | 0.586 | 0.901 |
| 9 | kNN cos_std | 0.570 | 0.576 | 0.450 | 0.767 |
| 9 | kNN cos_raw | 0.558 | 0.568 | 0.430 | 0.757 |
| 6 | linear (ridge) | 0.676 | 0.664 | 0.527 | 0.860 |
| 6 | MLP | 0.739 | 0.694 | 0.598 | 0.911 |
| 6 | kNN cos_std | 0.616 | 0.614 | 0.506 | 0.823 |
| 6 | kNN cos_raw | 0.612 | 0.614 | 0.491 | 0.809 |
| 3 | linear (ridge) | 0.675 | 0.662 | 0.521 | 0.856 |
| 3 | MLP | 0.733 | 0.692 | 0.600 | 0.910 |
| 3 | kNN cos_std | 0.613 | 0.597 | 0.500 | 0.824 |
| 3 | kNN cos_raw | 0.610 | 0.604 | 0.490 | 0.815 |
| in | linear (ridge) | 0.695 | 0.670 | 0.518 | 0.865 |
| in | MLP | 0.744 | 0.703 | 0.615 | 0.917 |
| in | kNN cos_std | 0.645 | 0.636 | 0.487 | 0.835 |
| in | kNN cos_raw | 0.648 | 0.637 | 0.490 | 0.836 |
| - | head | 0.693 | 0.681 | 0.564 | 0.890 |

## Gap ratio on the look-alike pairs: all (1506 pairs)

Median |pred NK - pred T| / median |measured NK - measured T|, [95% NK-cell bootstrap]; 1 = gap kept, 0 = gap gone. Last column: the same read-out on random NK x T pairs (median over 5 seeds), CD56 / CD94 / CD335 / CD3.

| depth | read-out | CD56 | CD94 | CD335 | CD3 | random pairs |
|---|---|---|---|---|---|---|
| z | linear (ridge) | 0.575 [0.520, 0.687] | 0.360 [0.299, 0.511] | 0.340 [0.315, 0.374] | 0.799 [0.718, 0.893] | 0.96 / 1.16 / 0.93 / 0.73 |
| z | MLP | 0.463 [0.422, 0.554] | 0.356 [0.298, 0.506] | 0.424 [0.395, 0.464] | 0.478 [0.402, 0.544] | 1.08 / 1.47 / 1.22 / 0.97 |
| z | kNN cos_std | 0.290 [0.256, 0.348] | 0.231 [0.187, 0.336] | 0.246 [0.219, 0.275] | 0.469 [0.397, 0.545] | 1.03 / 1.18 / 1.09 / 0.91 |
| z | kNN cos_raw | 0.345 [0.303, 0.418] | 0.227 [0.185, 0.330] | 0.227 [0.200, 0.257] | 0.526 [0.436, 0.611] | 1.01 / 1.17 / 1.05 / 0.87 |
| 12 | linear (ridge) | 0.560 [0.502, 0.662] | 0.337 [0.280, 0.481] | 0.333 [0.304, 0.365] | 0.806 [0.717, 0.896] | 0.96 / 1.15 / 0.93 / 0.73 |
| 12 | MLP | 0.471 [0.433, 0.556] | 0.339 [0.284, 0.479] | 0.422 [0.396, 0.462] | 0.491 [0.426, 0.558] | 1.08 / 1.44 / 1.25 / 0.94 |
| 12 | kNN cos_std | 0.291 [0.255, 0.357] | 0.224 [0.182, 0.318] | 0.244 [0.218, 0.275] | 0.470 [0.395, 0.565] | 1.03 / 1.18 / 1.09 / 0.91 |
| 12 | kNN cos_raw | 0.345 [0.303, 0.416] | 0.224 [0.184, 0.326] | 0.228 [0.201, 0.258] | 0.526 [0.437, 0.611] | 1.01 / 1.17 / 1.05 / 0.87 |
| 9 | linear (ridge) | 0.400 [0.356, 0.488] | 0.294 [0.244, 0.422] | 0.210 [0.190, 0.234] | 0.599 [0.536, 0.677] | 0.71 / 0.92 / 0.63 / 0.51 |
| 9 | MLP | 0.578 [0.509, 0.699] | 0.401 [0.334, 0.572] | 0.429 [0.382, 0.471] | 1.106 [0.986, 1.243] | 1.05 / 1.33 / 1.18 / 0.91 |
| 9 | kNN cos_std | 0.683 [0.595, 0.839] | 0.479 [0.391, 0.700] | 0.222 [0.197, 0.255] | 0.815 [0.697, 0.931] | 0.62 / 0.86 / 0.56 / 0.63 |
| 9 | kNN cos_raw | 0.623 [0.552, 0.766] | 0.430 [0.351, 0.620] | 0.177 [0.158, 0.207] | 0.810 [0.709, 0.916] | 0.54 / 0.70 / 0.43 / 0.54 |
| 6 | linear (ridge) | 0.493 [0.439, 0.597] | 0.356 [0.293, 0.507] | 0.256 [0.234, 0.282] | 0.646 [0.573, 0.720] | 0.78 / 0.98 / 0.71 / 0.56 |
| 6 | MLP | 0.570 [0.513, 0.676] | 0.418 [0.344, 0.594] | 0.436 [0.399, 0.480] | 1.094 [0.956, 1.232] | 1.04 / 1.30 / 1.17 / 0.92 |
| 6 | kNN cos_std | 0.576 [0.506, 0.709] | 0.397 [0.330, 0.581] | 0.236 [0.215, 0.265] | 0.834 [0.719, 0.965] | 0.76 / 1.04 / 0.78 / 0.67 |
| 6 | kNN cos_raw | 0.514 [0.451, 0.626] | 0.368 [0.308, 0.538] | 0.223 [0.203, 0.251] | 0.992 [0.862, 1.140] | 0.69 / 0.92 / 0.73 / 0.65 |
| 3 | linear (ridge) | 0.439 [0.393, 0.528] | 0.302 [0.249, 0.427] | 0.210 [0.191, 0.233] | 0.656 [0.579, 0.732] | 0.65 / 0.78 / 0.57 / 0.52 |
| 3 | MLP | 0.535 [0.477, 0.641] | 0.437 [0.365, 0.612] | 0.518 [0.472, 0.578] | 1.491 [1.317, 1.682] | 1.08 / 1.40 / 1.27 / 0.96 |
| 3 | kNN cos_std | 0.581 [0.508, 0.724] | 0.350 [0.288, 0.496] | 0.228 [0.201, 0.257] | 1.096 [0.952, 1.260] | 0.64 / 0.84 / 0.61 / 0.61 |
| 3 | kNN cos_raw | 0.509 [0.450, 0.629] | 0.286 [0.235, 0.410] | 0.212 [0.188, 0.237] | 0.977 [0.852, 1.123] | 0.61 / 0.76 / 0.61 / 0.61 |
| in | linear (ridge) | 0.436 [0.392, 0.525] | 0.293 [0.245, 0.414] | 0.216 [0.195, 0.241] | 0.653 [0.577, 0.729] | 0.64 / 0.75 / 0.55 / 0.51 |
| in | MLP | 0.516 [0.458, 0.606] | 0.413 [0.344, 0.580] | 0.548 [0.499, 0.602] | 1.290 [1.147, 1.442] | 1.12 / 1.38 / 1.25 / 0.97 |
| in | kNN cos_std | 0.425 [0.368, 0.515] | 0.383 [0.316, 0.554] | 0.261 [0.238, 0.293] | 1.005 [0.876, 1.159] | 0.66 / 0.96 / 0.71 / 0.60 |
| in | kNN cos_raw | 0.402 [0.359, 0.499] | 0.369 [0.305, 0.536] | 0.242 [0.218, 0.276] | 1.000 [0.867, 1.154] | 0.64 / 0.98 / 0.69 / 0.60 |
| - | head | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.115 [0.078, 0.147] | 0.154 [0.137, 0.179] | 1.18 / 1.52 / 1.43 / 0.77 |

MLP per-seed range of the gap ratio: z: CD56 0.428-0.528, CD94 0.365-0.374, CD335 0.418-0.435, CD3 0.426-0.480; 12: CD56 0.453-0.520, CD94 0.268-0.411, CD335 0.411-0.434, CD3 0.246-0.712; 9: CD56 0.578-0.657, CD94 0.405-0.450, CD335 0.377-0.477, CD3 0.744-1.389; 6: CD56 0.604-0.644, CD94 0.433-0.466, CD335 0.440-0.455, CD3 1.034-1.056; 3: CD56 0.502-0.648, CD94 0.426-0.474, CD335 0.490-0.599, CD3 1.204-1.691; in: CD56 0.467-0.597, CD94 0.390-0.450, CD335 0.520-0.555, CD3 1.158-1.479.

### Paired differences in gap ratio (all); interval above 0 = the first read-out keeps more of the gap

| depth | comparison | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|---|
| z | MLP - linear | -0.112 [-0.175, -0.057] | -0.004 [-0.041, +0.031] | +0.084 [+0.061, +0.108] | -0.322 [-0.440, -0.222] |
| z | MLP - head | +0.463 [+0.422, +0.554] | +0.356 [+0.298, +0.506] | +0.309 [+0.277, +0.350] | +0.323 [+0.256, +0.377] |
| z | kNN cos_std - linear | -0.285 [-0.358, -0.236] | -0.129 [-0.188, -0.097] | -0.094 [-0.120, -0.067] | -0.331 [-0.434, -0.234] |
| z | kNN cos_std - head | +0.290 [+0.256, +0.348] | +0.231 [+0.187, +0.336] | +0.131 [+0.099, +0.169] | +0.314 [+0.253, +0.382] |
| z | kNN cos_raw - linear | -0.230 [-0.303, -0.181] | -0.133 [-0.189, -0.099] | -0.113 [-0.141, -0.087] | -0.274 [-0.386, -0.182] |
| z | kNN cos_raw - head | +0.345 [+0.303, +0.418] | +0.227 [+0.185, +0.330] | +0.112 [+0.080, +0.152] | +0.371 [+0.287, +0.452] |
| 12 | MLP - linear | -0.089 [-0.144, -0.041] | +0.001 [-0.030, +0.028] | +0.089 [+0.068, +0.114] | -0.316 [-0.408, -0.209] |
| 12 | MLP - head | +0.471 [+0.433, +0.556] | +0.339 [+0.284, +0.479] | +0.307 [+0.278, +0.349] | +0.337 [+0.281, +0.392] |
| 12 | kNN cos_std - linear | -0.269 [-0.336, -0.211] | -0.114 [-0.169, -0.086] | -0.089 [-0.113, -0.057] | -0.337 [-0.433, -0.220] |
| 12 | kNN cos_std - head | +0.291 [+0.255, +0.357] | +0.224 [+0.182, +0.318] | +0.129 [+0.099, +0.169] | +0.316 [+0.250, +0.401] |
| 12 | kNN cos_raw - linear | -0.215 [-0.283, -0.161] | -0.113 [-0.166, -0.082] | -0.105 [-0.133, -0.078] | -0.281 [-0.383, -0.178] |
| 12 | kNN cos_raw - head | +0.345 [+0.303, +0.416] | +0.224 [+0.184, +0.326] | +0.113 [+0.080, +0.153] | +0.371 [+0.287, +0.452] |
| 9 | MLP - linear | +0.178 [+0.127, +0.236] | +0.107 [+0.072, +0.158] | +0.218 [+0.178, +0.250] | +0.507 [+0.412, +0.600] |
| 9 | MLP - head | +0.578 [+0.509, +0.699] | +0.401 [+0.334, +0.572] | +0.314 [+0.266, +0.357] | +0.951 [+0.833, +1.082] |
| 9 | kNN cos_std - linear | +0.283 [+0.211, +0.378] | +0.185 [+0.128, +0.280] | +0.012 [-0.011, +0.042] | +0.216 [+0.117, +0.297] |
| 9 | kNN cos_std - head | +0.683 [+0.595, +0.839] | +0.479 [+0.391, +0.700] | +0.107 [+0.071, +0.158] | +0.661 [+0.550, +0.767] |
| 9 | kNN cos_raw - linear | +0.223 [+0.162, +0.295] | +0.136 [+0.093, +0.209] | -0.033 [-0.055, -0.007] | +0.211 [+0.124, +0.286] |
| 9 | kNN cos_raw - head | +0.623 [+0.552, +0.766] | +0.430 [+0.351, +0.620] | +0.062 [+0.028, +0.106] | +0.655 [+0.556, +0.759] |
| 6 | MLP - linear | +0.077 [+0.029, +0.126] | +0.062 [+0.023, +0.104] | +0.180 [+0.151, +0.211] | +0.448 [+0.334, +0.564] |
| 6 | MLP - head | +0.570 [+0.513, +0.676] | +0.418 [+0.344, +0.594] | +0.321 [+0.282, +0.370] | +0.939 [+0.803, +1.073] |
| 6 | kNN cos_std - linear | +0.083 [+0.024, +0.148] | +0.040 [+0.003, +0.098] | -0.020 [-0.041, +0.005] | +0.187 [+0.106, +0.284] |
| 6 | kNN cos_std - head | +0.576 [+0.506, +0.709] | +0.397 [+0.330, +0.581] | +0.121 [+0.087, +0.168] | +0.679 [+0.565, +0.806] |
| 6 | kNN cos_raw - linear | +0.021 [-0.037, +0.076] | +0.012 [-0.025, +0.057] | -0.034 [-0.054, -0.007] | +0.345 [+0.255, +0.452] |
| 6 | kNN cos_raw - head | +0.514 [+0.451, +0.626] | +0.368 [+0.308, +0.538] | +0.107 [+0.074, +0.154] | +0.837 [+0.704, +0.981] |
| 3 | MLP - linear | +0.096 [+0.046, +0.151] | +0.134 [+0.094, +0.199] | +0.309 [+0.270, +0.355] | +0.835 [+0.682, +0.985] |
| 3 | MLP - head | +0.535 [+0.477, +0.641] | +0.437 [+0.365, +0.612] | +0.404 [+0.355, +0.463] | +1.337 [+1.162, +1.529] |
| 3 | kNN cos_std - linear | +0.142 [+0.089, +0.219] | +0.048 [+0.006, +0.097] | +0.018 [-0.006, +0.042] | +0.440 [+0.338, +0.561] |
| 3 | kNN cos_std - head | +0.581 [+0.508, +0.724] | +0.350 [+0.288, +0.496] | +0.112 [+0.074, +0.157] | +0.942 [+0.798, +1.104] |
| 3 | kNN cos_raw - linear | +0.070 [+0.025, +0.125] | -0.017 [-0.053, +0.017] | +0.002 [-0.021, +0.022] | +0.320 [+0.239, +0.427] |
| 3 | kNN cos_raw - head | +0.509 [+0.450, +0.629] | +0.286 [+0.235, +0.410] | +0.097 [+0.061, +0.138] | +0.822 [+0.702, +0.966] |
| in | MLP - linear | +0.080 [+0.022, +0.127] | +0.120 [+0.081, +0.178] | +0.331 [+0.292, +0.371] | +0.637 [+0.527, +0.748] |
| in | MLP - head | +0.516 [+0.458, +0.606] | +0.413 [+0.344, +0.580] | +0.433 [+0.386, +0.491] | +1.135 [+0.993, +1.282] |
| in | kNN cos_std - linear | -0.011 [-0.062, +0.036] | +0.090 [+0.049, +0.147] | +0.044 [+0.021, +0.074] | +0.352 [+0.249, +0.481] |
| in | kNN cos_std - head | +0.425 [+0.368, +0.515] | +0.383 [+0.316, +0.554] | +0.146 [+0.112, +0.192] | +0.851 [+0.719, +1.002] |
| in | kNN cos_raw - linear | -0.034 [-0.076, +0.017] | +0.076 [+0.036, +0.132] | +0.026 [+0.004, +0.051] | +0.347 [+0.242, +0.470] |
| in | kNN cos_raw - head | +0.402 [+0.359, +0.499] | +0.369 [+0.305, +0.536] | +0.127 [+0.094, +0.170] | +0.846 [+0.713, +0.996] |

## Gap ratio on the look-alike pairs: no_gdT158 (453 pairs)

Median |pred NK - pred T| / median |measured NK - measured T|, [95% NK-cell bootstrap]; 1 = gap kept, 0 = gap gone. Last column: the same read-out on random NK x T pairs (median over 5 seeds), CD56 / CD94 / CD335 / CD3.

| depth | read-out | CD56 | CD94 | CD335 | CD3 | random pairs |
|---|---|---|---|---|---|---|
| z | linear (ridge) | 0.279 [0.239, 0.316] | 0.397 [0.327, 0.498] | 0.404 [0.343, 0.460] | 0.276 [0.231, 0.304] | 0.97 / 1.17 / 0.95 / 0.73 |
| z | MLP | 0.471 [0.395, 0.554] | 0.589 [0.492, 0.732] | 0.740 [0.623, 0.830] | 0.609 [0.495, 0.673] | 1.07 / 1.46 / 1.24 / 0.96 |
| z | kNN cos_std | 0.215 [0.158, 0.267] | 0.272 [0.208, 0.365] | 0.343 [0.268, 0.436] | 0.280 [0.212, 0.344] | 1.03 / 1.19 / 1.11 / 0.91 |
| z | kNN cos_raw | 0.223 [0.175, 0.263] | 0.287 [0.216, 0.387] | 0.308 [0.233, 0.386] | 0.286 [0.221, 0.345] | 1.01 / 1.18 / 1.06 / 0.88 |
| 12 | linear (ridge) | 0.265 [0.228, 0.301] | 0.365 [0.307, 0.468] | 0.390 [0.336, 0.446] | 0.264 [0.224, 0.293] | 0.96 / 1.16 / 0.95 / 0.73 |
| 12 | MLP | 0.482 [0.405, 0.546] | 0.589 [0.481, 0.744] | 0.753 [0.615, 0.840] | 0.547 [0.465, 0.625] | 1.07 / 1.44 / 1.26 / 0.94 |
| 12 | kNN cos_std | 0.207 [0.162, 0.266] | 0.274 [0.207, 0.367] | 0.335 [0.262, 0.430] | 0.262 [0.202, 0.332] | 1.03 / 1.19 / 1.11 / 0.91 |
| 12 | kNN cos_raw | 0.224 [0.175, 0.263] | 0.289 [0.216, 0.386] | 0.308 [0.233, 0.387] | 0.286 [0.221, 0.345] | 1.01 / 1.18 / 1.06 / 0.88 |
| 9 | linear (ridge) | 0.155 [0.128, 0.191] | 0.231 [0.190, 0.305] | 0.217 [0.176, 0.254] | 0.188 [0.160, 0.218] | 0.72 / 0.93 / 0.65 / 0.52 |
| 9 | MLP | 0.339 [0.271, 0.402] | 0.466 [0.372, 0.595] | 0.540 [0.453, 0.645] | 0.415 [0.343, 0.503] | 1.05 / 1.32 / 1.20 / 0.92 |
| 9 | kNN cos_std | 0.184 [0.155, 0.223] | 0.310 [0.250, 0.408] | 0.161 [0.130, 0.194] | 0.152 [0.124, 0.212] | 0.63 / 0.86 / 0.57 / 0.64 |
| 9 | kNN cos_raw | 0.175 [0.148, 0.208] | 0.301 [0.249, 0.384] | 0.145 [0.115, 0.170] | 0.188 [0.149, 0.217] | 0.55 / 0.71 / 0.44 / 0.56 |
| 6 | linear (ridge) | 0.179 [0.154, 0.218] | 0.296 [0.239, 0.372] | 0.255 [0.218, 0.291] | 0.213 [0.181, 0.240] | 0.79 / 1.00 / 0.73 / 0.57 |
| 6 | MLP | 0.366 [0.289, 0.451] | 0.475 [0.384, 0.608] | 0.603 [0.505, 0.696] | 0.484 [0.389, 0.563] | 1.03 / 1.29 / 1.18 / 0.92 |
| 6 | kNN cos_std | 0.163 [0.139, 0.188] | 0.290 [0.218, 0.377] | 0.182 [0.146, 0.214] | 0.162 [0.129, 0.194] | 0.78 / 1.05 / 0.80 / 0.69 |
| 6 | kNN cos_raw | 0.154 [0.125, 0.178] | 0.295 [0.238, 0.382] | 0.192 [0.157, 0.228] | 0.200 [0.164, 0.240] | 0.70 / 0.93 / 0.75 / 0.66 |
| 3 | linear (ridge) | 0.185 [0.153, 0.219] | 0.278 [0.225, 0.354] | 0.241 [0.198, 0.280] | 0.203 [0.179, 0.230] | 0.65 / 0.80 / 0.59 / 0.53 |
| 3 | MLP | 0.364 [0.290, 0.464] | 0.551 [0.447, 0.711] | 0.751 [0.597, 0.875] | 0.501 [0.426, 0.569] | 1.07 / 1.40 / 1.28 / 0.96 |
| 3 | kNN cos_std | 0.142 [0.113, 0.178] | 0.253 [0.201, 0.322] | 0.175 [0.135, 0.214] | 0.210 [0.169, 0.258] | 0.65 / 0.85 / 0.62 / 0.62 |
| 3 | kNN cos_raw | 0.136 [0.110, 0.168] | 0.208 [0.166, 0.271] | 0.178 [0.141, 0.218] | 0.186 [0.161, 0.218] | 0.62 / 0.77 / 0.63 / 0.62 |
| in | linear (ridge) | 0.218 [0.182, 0.254] | 0.302 [0.254, 0.386] | 0.293 [0.241, 0.331] | 0.207 [0.188, 0.236] | 0.64 / 0.76 / 0.57 / 0.52 |
| in | MLP | 0.433 [0.356, 0.523] | 0.575 [0.486, 0.709] | 0.796 [0.680, 0.913] | 0.549 [0.484, 0.616] | 1.12 / 1.38 / 1.26 / 0.97 |
| in | kNN cos_std | 0.143 [0.110, 0.172] | 0.305 [0.248, 0.396] | 0.280 [0.228, 0.341] | 0.228 [0.198, 0.260] | 0.67 / 0.97 / 0.73 / 0.62 |
| in | kNN cos_raw | 0.146 [0.116, 0.176] | 0.295 [0.239, 0.384] | 0.270 [0.217, 0.334] | 0.241 [0.204, 0.277] | 0.66 / 0.99 / 0.71 / 0.62 |
| - | head | 0.238 [0.167, 0.312] | 0.311 [0.199, 0.427] | 0.880 [0.695, 1.009] | 0.271 [0.199, 0.364] | 1.18 / 1.53 / 1.44 / 0.77 |

MLP per-seed range of the gap ratio: z: CD56 0.423-0.492, CD94 0.554-0.625, CD335 0.723-0.746, CD3 0.561-0.641; 12: CD56 0.418-0.505, CD94 0.533-0.679, CD335 0.728-0.762, CD3 0.514-0.586; 9: CD56 0.300-0.351, CD94 0.427-0.479, CD335 0.516-0.605, CD3 0.346-0.457; 6: CD56 0.349-0.413, CD94 0.440-0.562, CD335 0.590-0.670, CD3 0.461-0.501; 3: CD56 0.340-0.387, CD94 0.525-0.615, CD335 0.657-0.764, CD3 0.509-0.517; in: CD56 0.398-0.461, CD94 0.545-0.589, CD335 0.759-0.848, CD3 0.540-0.572.

### Paired differences in gap ratio (no_gdT158); interval above 0 = the first read-out keeps more of the gap

| depth | comparison | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|---|
| z | MLP - linear | +0.192 [+0.122, +0.266] | +0.192 [+0.100, +0.284] | +0.336 [+0.241, +0.406] | +0.333 [+0.237, +0.389] |
| z | MLP - head | +0.233 [+0.182, +0.298] | +0.278 [+0.202, +0.382] | -0.140 [-0.225, -0.026] | +0.338 [+0.242, +0.393] |
| z | kNN cos_std - linear | -0.064 [-0.115, -0.017] | -0.124 [-0.186, -0.072] | -0.061 [-0.140, +0.042] | +0.003 [-0.061, +0.071] |
| z | kNN cos_std - head | -0.023 [-0.090, +0.034] | -0.039 [-0.118, +0.081] | -0.537 [-0.641, -0.382] | +0.008 [-0.079, +0.070] |
| z | kNN cos_raw - linear | -0.056 [-0.099, -0.015] | -0.109 [-0.171, -0.052] | -0.096 [-0.180, -0.011] | +0.010 [-0.054, +0.078] |
| z | kNN cos_raw - head | -0.015 [-0.084, +0.043] | -0.024 [-0.105, +0.099] | -0.572 [-0.677, -0.419] | +0.015 [-0.089, +0.089] |
| 12 | MLP - linear | +0.217 [+0.150, +0.273] | +0.224 [+0.120, +0.321] | +0.363 [+0.244, +0.437] | +0.283 [+0.219, +0.354] |
| 12 | MLP - head | +0.244 [+0.186, +0.287] | +0.278 [+0.204, +0.382] | -0.127 [-0.205, -0.034] | +0.276 [+0.196, +0.333] |
| 12 | kNN cos_std - linear | -0.057 [-0.104, -0.005] | -0.091 [-0.163, -0.035] | -0.055 [-0.128, +0.044] | -0.002 [-0.059, +0.074] |
| 12 | kNN cos_std - head | -0.031 [-0.089, +0.030] | -0.037 [-0.118, +0.084] | -0.545 [-0.649, -0.391] | -0.009 [-0.092, +0.067] |
| 12 | kNN cos_raw - linear | -0.040 [-0.084, -0.005] | -0.076 [-0.144, -0.028] | -0.083 [-0.158, -0.003] | +0.022 [-0.043, +0.087] |
| 12 | kNN cos_raw - head | -0.013 [-0.084, +0.043] | -0.022 [-0.105, +0.100] | -0.572 [-0.676, -0.418] | +0.015 [-0.089, +0.088] |
| 9 | MLP - linear | +0.184 [+0.130, +0.235] | +0.236 [+0.150, +0.318] | +0.323 [+0.259, +0.411] | +0.227 [+0.173, +0.297] |
| 9 | MLP - head | +0.101 [+0.042, +0.159] | +0.155 [+0.062, +0.268] | -0.340 [-0.418, -0.186] | +0.144 [+0.053, +0.230] |
| 9 | kNN cos_std - linear | +0.029 [-0.009, +0.073] | +0.080 [+0.015, +0.154] | -0.056 [-0.095, -0.010] | -0.036 [-0.067, +0.015] |
| 9 | kNN cos_std - head | -0.054 [-0.111, +0.011] | -0.001 [-0.086, +0.121] | -0.719 [-0.838, -0.540] | -0.120 [-0.217, -0.034] |
| 9 | kNN cos_raw - linear | +0.020 [-0.016, +0.061] | +0.070 [+0.018, +0.124] | -0.071 [-0.113, -0.038] | +0.000 [-0.036, +0.031] |
| 9 | kNN cos_raw - head | -0.063 [-0.127, +0.003] | -0.011 [-0.090, +0.113] | -0.734 [-0.859, -0.556] | -0.083 [-0.178, -0.011] |
| 6 | MLP - linear | +0.187 [+0.115, +0.254] | +0.180 [+0.101, +0.266] | +0.348 [+0.264, +0.429] | +0.271 [+0.200, +0.343] |
| 6 | MLP - head | +0.128 [+0.069, +0.194] | +0.164 [+0.081, +0.281] | -0.277 [-0.360, -0.140] | +0.213 [+0.108, +0.292] |
| 6 | kNN cos_std - linear | -0.016 [-0.056, +0.009] | -0.006 [-0.073, +0.052] | -0.073 [-0.116, -0.038] | -0.051 [-0.086, -0.018] |
| 6 | kNN cos_std - head | -0.075 [-0.140, -0.012] | -0.021 [-0.124, +0.096] | -0.698 [-0.820, -0.523] | -0.110 [-0.207, -0.042] |
| 6 | kNN cos_raw - linear | -0.025 [-0.066, -0.000] | -0.001 [-0.045, +0.056] | -0.063 [-0.102, -0.028] | -0.013 [-0.052, +0.030] |
| 6 | kNN cos_raw - head | -0.084 [-0.147, -0.018] | -0.016 [-0.094, +0.104] | -0.688 [-0.808, -0.510] | -0.071 [-0.169, +0.005] |
| 3 | MLP - linear | +0.179 [+0.122, +0.261] | +0.273 [+0.190, +0.394] | +0.510 [+0.367, +0.637] | +0.299 [+0.238, +0.355] |
| 3 | MLP - head | +0.126 [+0.070, +0.200] | +0.239 [+0.148, +0.371] | -0.129 [-0.237, +0.011] | +0.230 [+0.151, +0.299] |
| 3 | kNN cos_std - linear | -0.043 [-0.077, -0.005] | -0.025 [-0.076, +0.022] | -0.067 [-0.105, -0.027] | +0.008 [-0.034, +0.053] |
| 3 | kNN cos_std - head | -0.096 [-0.170, -0.020] | -0.058 [-0.146, +0.060] | -0.705 [-0.831, -0.527] | -0.061 [-0.169, +0.023] |
| 3 | kNN cos_raw - linear | -0.049 [-0.083, -0.016] | -0.070 [-0.118, -0.020] | -0.063 [-0.101, -0.022] | -0.017 [-0.048, +0.017] |
| 3 | kNN cos_raw - head | -0.102 [-0.171, -0.032] | -0.103 [-0.194, +0.015] | -0.702 [-0.827, -0.527] | -0.086 [-0.186, -0.005] |
| in | MLP - linear | +0.214 [+0.151, +0.285] | +0.273 [+0.198, +0.356] | +0.503 [+0.416, +0.619] | +0.343 [+0.286, +0.400] |
| in | MLP - head | +0.195 [+0.139, +0.260] | +0.264 [+0.184, +0.384] | -0.084 [-0.187, +0.058] | +0.278 [+0.194, +0.339] |
| in | kNN cos_std - linear | -0.076 [-0.115, -0.041] | +0.003 [-0.055, +0.064] | -0.013 [-0.056, +0.055] | +0.021 [-0.011, +0.050] |
| in | kNN cos_std - head | -0.095 [-0.168, -0.028] | -0.006 [-0.091, +0.119] | -0.600 [-0.716, -0.422] | -0.043 [-0.128, +0.022] |
| in | kNN cos_raw - linear | -0.072 [-0.109, -0.036] | -0.006 [-0.073, +0.057] | -0.024 [-0.059, +0.044] | +0.034 [-0.004, +0.066] |
| in | kNN cos_raw - head | -0.091 [-0.150, -0.030] | -0.016 [-0.104, +0.106] | -0.610 [-0.722, -0.437] | -0.030 [-0.122, +0.037] |

## Look-alike gap ratio divided by the same read-out's random-pair ratio

A read-out that shrinks every prediction toward the mean lowers both ratios; this quotient asks how much MORE the look-alike gap is compressed than a random NK x T gap by the same read-out (1 = no look-alike-specific loss). Point values only (derived from the two tables above).

| pair set | depth | read-out | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|---|---|
| all | z | linear (ridge) | 0.599 | 0.309 | 0.365 | 1.098 |
| all | z | MLP | 0.428 | 0.242 | 0.346 | 0.493 |
| all | z | kNN cos_std | 0.282 | 0.196 | 0.225 | 0.516 |
| all | z | kNN cos_raw | 0.340 | 0.194 | 0.217 | 0.604 |
| all | 12 | linear (ridge) | 0.586 | 0.293 | 0.358 | 1.107 |
| all | 12 | MLP | 0.437 | 0.235 | 0.337 | 0.522 |
| all | 12 | kNN cos_std | 0.283 | 0.189 | 0.223 | 0.518 |
| all | 12 | kNN cos_raw | 0.340 | 0.192 | 0.217 | 0.604 |
| all | 9 | linear (ridge) | 0.565 | 0.320 | 0.333 | 1.177 |
| all | 9 | MLP | 0.551 | 0.302 | 0.364 | 1.211 |
| all | 9 | kNN cos_std | 1.097 | 0.559 | 0.400 | 1.304 |
| all | 9 | kNN cos_raw | 1.148 | 0.610 | 0.413 | 1.503 |
| all | 6 | linear (ridge) | 0.632 | 0.365 | 0.361 | 1.148 |
| all | 6 | MLP | 0.549 | 0.323 | 0.375 | 1.189 |
| all | 6 | kNN cos_std | 0.760 | 0.383 | 0.304 | 1.240 |
| all | 6 | kNN cos_raw | 0.751 | 0.402 | 0.305 | 1.532 |
| all | 3 | linear (ridge) | 0.679 | 0.386 | 0.366 | 1.268 |
| all | 3 | MLP | 0.497 | 0.311 | 0.409 | 1.552 |
| all | 3 | kNN cos_std | 0.914 | 0.418 | 0.373 | 1.790 |
| all | 3 | kNN cos_raw | 0.832 | 0.375 | 0.346 | 1.610 |
| all | in | linear (ridge) | 0.684 | 0.390 | 0.392 | 1.290 |
| all | in | MLP | 0.460 | 0.299 | 0.440 | 1.330 |
| all | in | kNN cos_std | 0.648 | 0.397 | 0.369 | 1.673 |
| all | in | kNN cos_raw | 0.625 | 0.378 | 0.350 | 1.656 |
| all | - | head | 0.000 | 0.000 | 0.081 | 0.201 |
| no_gdT158 | z | linear (ridge) | 0.289 | 0.339 | 0.427 | 0.378 |
| no_gdT158 | z | MLP | 0.439 | 0.403 | 0.599 | 0.631 |
| no_gdT158 | z | kNN cos_std | 0.209 | 0.229 | 0.310 | 0.306 |
| no_gdT158 | z | kNN cos_raw | 0.220 | 0.243 | 0.290 | 0.326 |
| no_gdT158 | 12 | linear (ridge) | 0.275 | 0.315 | 0.412 | 0.360 |
| no_gdT158 | 12 | MLP | 0.451 | 0.410 | 0.595 | 0.584 |
| no_gdT158 | 12 | kNN cos_std | 0.201 | 0.230 | 0.303 | 0.288 |
| no_gdT158 | 12 | kNN cos_raw | 0.221 | 0.245 | 0.290 | 0.327 |
| no_gdT158 | 9 | linear (ridge) | 0.216 | 0.249 | 0.334 | 0.360 |
| no_gdT158 | 9 | MLP | 0.323 | 0.354 | 0.451 | 0.454 |
| no_gdT158 | 9 | kNN cos_std | 0.291 | 0.360 | 0.281 | 0.235 |
| no_gdT158 | 9 | kNN cos_raw | 0.318 | 0.421 | 0.330 | 0.336 |
| no_gdT158 | 6 | linear (ridge) | 0.226 | 0.297 | 0.349 | 0.370 |
| no_gdT158 | 6 | MLP | 0.355 | 0.368 | 0.510 | 0.525 |
| no_gdT158 | 6 | kNN cos_std | 0.210 | 0.276 | 0.228 | 0.236 |
| no_gdT158 | 6 | kNN cos_raw | 0.219 | 0.317 | 0.257 | 0.303 |
| no_gdT158 | 3 | linear (ridge) | 0.283 | 0.348 | 0.410 | 0.384 |
| no_gdT158 | 3 | MLP | 0.340 | 0.394 | 0.586 | 0.522 |
| no_gdT158 | 3 | kNN cos_std | 0.220 | 0.299 | 0.281 | 0.338 |
| no_gdT158 | 3 | kNN cos_raw | 0.219 | 0.270 | 0.283 | 0.300 |
| no_gdT158 | in | linear (ridge) | 0.339 | 0.396 | 0.516 | 0.402 |
| no_gdT158 | in | MLP | 0.387 | 0.418 | 0.633 | 0.566 |
| no_gdT158 | in | kNN cos_std | 0.213 | 0.314 | 0.386 | 0.371 |
| no_gdT158 | in | kNN cos_raw | 0.223 | 0.299 | 0.379 | 0.390 |
| no_gdT158 | - | head | 0.202 | 0.204 | 0.612 | 0.353 |

## kNN construction bias: shared training neighbours

Mean share of the k nearest training cells that the two cells of a pair have in common (k = the protein's chosen k; CD56 shown, other proteins in the JSON). High for look-alike pairs at z by construction.

| depth | metric | k (CD56) | look-alike all | look-alike no_gdT158 | random NK x T |
|---|---|---|---|---|---|
| z | cos_std | 100 | 0.443 | 0.384 | 0.007 |
| z | cos_raw | 100 | 0.439 | 0.379 | 0.007 |
| 12 | cos_std | 100 | 0.441 | 0.386 | 0.007 |
| 12 | cos_raw | 100 | 0.439 | 0.379 | 0.007 |
| 9 | cos_std | 20 | 0.125 | 0.139 | 0.006 |
| 9 | cos_raw | 20 | 0.111 | 0.117 | 0.005 |
| 6 | cos_std | 50 | 0.143 | 0.153 | 0.011 |
| 6 | cos_raw | 50 | 0.163 | 0.158 | 0.016 |
| 3 | cos_std | 50 | 0.116 | 0.139 | 0.008 |
| 3 | cos_raw | 50 | 0.142 | 0.145 | 0.013 |
| in | cos_std | 50 | 0.096 | 0.102 | 0.010 |
| in | cos_raw | 50 | 0.084 | 0.087 | 0.010 |

## Reading

**Short answer.** On the clean look-alike pairs (no_gdT158, 453 pairs), a non-linear read-out (the MLP)
of TEDDY's pooled layers keeps clearly more of the NK-T protein gap than the linear probe does, at every
depth tested and for all four proteins. For CD56, CD94 and CD3 it also keeps more than our head. For
CD335 it does not beat the head, which already keeps 0.880. By the rule for this step, the difference is
in the pooled (gene-mean) representation in a form a linear read-out misses, and a better read-out gets
part of it back. It does not get all of it back, and the 'all' pair set does not support the same
statement.

At z, no_gdT158 (gap ratio [95% CI]):

| protein | MLP | linear | head | MLP - linear | MLP - head |
|---|---|---|---|---|---|
| CD56 | 0.471 [0.395, 0.554] | 0.279 [0.239, 0.316] | 0.238 [0.167, 0.312] | +0.192 [+0.122, +0.266] | +0.233 [+0.182, +0.298] |
| CD94 | 0.589 [0.492, 0.732] | 0.397 [0.327, 0.498] | 0.311 [0.199, 0.427] | +0.192 [+0.100, +0.284] | +0.278 [+0.202, +0.382] |
| CD335 | 0.740 [0.623, 0.830] | 0.404 [0.343, 0.460] | 0.880 [0.695, 1.009] | +0.336 [+0.241, +0.406] | -0.140 [-0.225, -0.026] |
| CD3 | 0.609 [0.495, 0.673] | 0.276 [0.231, 0.304] | 0.271 [0.199, 0.364] | +0.333 [+0.237, +0.389] | +0.338 [+0.242, +0.393] |

The MLP - linear interval is above 0 for all four proteins at all six depths (in, 3, 6, 9, 12, z). The
MLP - head interval is above 0 for CD56, CD94 and CD3 at all six depths. The MLP's site4 Pearson is also
higher than the ridge's and the head's at every depth (z: 0.779 / 0.712 / 0.613 / 0.932, against ridge
0.718 / 0.687 / 0.582 / 0.896 and head 0.693 / 0.681 / 0.564 / 0.890).

**What limits this.**
1. *Part of the gain is not specific to look-alikes.* The MLP shrinks its predictions less than the ridge,
   so its random NK x T ratios are higher too (z: 1.07 / 1.46 / 1.24 / 0.96, against ridge 0.97 / 1.17 /
   0.95 / 0.73). Divided by the random-pair ratio, the MLP still keeps more at z (0.44 / 0.40 / 0.60 /
   0.63, against ridge 0.29 / 0.34 / 0.43 / 0.38 and head 0.20 / 0.20 / 0.61 / 0.35). The gain is
   smallest for CD94. These quotients are point values with no interval.
2. *About half of the gap is still missing.* At z the MLP keeps 0.47 to 0.74 of the measured median gap,
   and its look-alike gap is still compressed relative to random pairs (quotient 0.40 to 0.63, below 1).
   A ratio of 1 is not the reachable ceiling, because the measured gap also contains ADT noise that no
   RNA-based read-out can predict. This run has no estimate of that ceiling.
3. *Depth does not add to it.* The MLP on the input embedding (the gene-mean of token + position
   embeddings, that is, which genes are among the cell's top tokens) keeps about as much as at z
   (no_gdT158: 0.433 / 0.575 / 0.796 / 0.549, against z 0.471 / 0.589 / 0.740 / 0.609). Layers 6 and 9
   keep less than both the input and z for all four proteins (point values; the intervals overlap). So what
   the MLP reads is mostly present in the cell's gene set already. This run shows no NK-T information that
   TEDDY's layers add on top of it for this read-out.
4. *The 'all' pair set is mixed.* 1,053 of its 1,506 pairs involve a gdT CD158b+ cell, and its measured
   CD3 gap is small (median 0.137). At z and layer 12 the MLP keeps *less* CD56 and CD3 gap than the ridge
   (z CD3: 0.478 against 0.799, MLP - linear -0.322 [-0.440, -0.222]) and more CD335. At the input and at
   layers 3 to 9 the MLP overshoots the measured CD3 gap (ratio 1.09 to 1.49). The MLP beats the head for
   all four proteins at every depth there (head 0.000 / 0.000 / 0.115 / 0.154). But "non-linear keeps
   more than linear" does not hold for this set.
5. *kNN does not recover the gap at any depth, and at z it is biased by construction.* The pairs were
   selected as neighbours at z, so the two cells of a look-alike pair share on average 0.38 to 0.44 of
   their 100 nearest training cells at z and layer 12 (k = 100, CD56's k there; random NK x T pairs: 0.007). They get near-identical
   kNN predictions, and kNN keeps less gap than the ridge there (no_gdT158, z, cos_raw: 0.223 / 0.287 /
   0.308 / 0.286). At the earlier depths the pairs were not selected and share only 0.08 to 0.16 of their
   neighbours (at CD56's k). On no_gdT158, kNN is still no better than the ridge at those depths: kNN - linear lies below
   0 or includes 0, except CD94 at layer 9 (cos_std +0.080 [+0.015, +0.154], cos_raw +0.070 [+0.018, +0.124]). kNN is a local average, so
   it cannot separate two cells whose neighbourhoods overlap. Its failure says little beyond that.
6. *For our head.* The head is itself an MLP on z, yet it keeps much less CD56, CD94 and CD3 gap than the
   MLP probe on the same z. The two differ in objective and target: the probe uses squared error on four
   p95-scaled proteins, and the head uses an NB likelihood on 134 protein counts with a constant size
   factor. This run does not test which of those differences costs the gap.
7. *Scope.* One held-out site, one pair construction (k = 10), one MLP architecture (2 x 256, early
   stopping at the best of epochs 6 to 18, patience 8), and 3 seeds with visible spread (z CD56 on
   no_gdT158: 0.423 to 0.492). There was no hyperparameter search beyond early stopping and the kNN k.

**Bottom line.** The pooled representation is not the whole bottleneck. On the clean pairs, a better
read-out of the same z keeps about 0.19 to 0.34 more of the measured gap than the ridge, and about 0.23 to
0.34 more than our head for CD56, CD94 and CD3. About half of the look-alike gap stays missing for every
read-out and depth tested, and the MLP on the input embedding does as well as on z. So it is still open
whether the per-token states carry the rest. `scripts/mode_a_token_probes.py` is written for that
question and has not been run.

## Token-level follow-up (written, not run)

`scripts/mode_a_token_probes.py` re-embeds 1,473 site4 cells with the official preprocessing: the 483
cells of the no_gdT158 pairs, plus 990 other site4 NK / T cells stratified to the pairs' token counts.
For the input, layer 6 and layer 12 it keeps the states of 13 gene tokens: NCAM1, KLRD1, NCR1, CD3E,
CD3D, CD3G, CD5, CD6, CD28, KLRF1, FCGR3A, KLRC1 and SH2D1B. The TCR constant genes are not in the
processed genes or the TEDDY vocab. It also keeps the gene-mean and a CLS-free attention-pooled summary.
It then fits ridge, MLP and kNN within site4 (trained on the random cells, scored on the pair cells),
against the gene-mean of the same forward pass. Command (repo root):

    PY=<venv>/bin/python
    $PY scripts/mode_a_token_probes.py --stage all \
        --processed data/processed/cite --embed-dir data/processed/cite_official \
        --ckpt ../teddy_mwe/ckpt/teddy_g_70M --medians data/reference/teddy_gene_medians.json \
        --linear-dir outputs/mode_a_official --out-dir outputs/mode_a_official/token_probes \
        --device auto --autocast fp16 --batch-size 4 --threads 4

GPU time estimate (MPS, as the official run): about 3 to 4 min for the forward pass. This is scaled from
the official run's measured 0.118 s/cell; the selected cells average 1,475 tokens against 1,266 overall.
Add the slower attention-weight path at 3 layers and batch 4, plus model load and tokenisation, and
expect under 10 min wall time. The CPU fit stage took 20 s in a smoke test on placeholder states. A
4-cell CPU smoke test (fp32, 4 threads, 3.1 s for cells of 1,696 to 2,048 tokens) matched the module's
own layer forward exactly (max |diff| 0.0). Its gene-means of layers 6 and 12 matched
z_rna_layer_means.npy to within 0.00024 and 0.0039. On CPU only, the full embed would take roughly 11 to
20 min, extrapolated from those 4 long cells.
