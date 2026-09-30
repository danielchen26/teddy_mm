# Mode A: layer-wise NK-T probes

Written by `scripts/mode_a_layer_probes.py` (mode_a_layer_probes v1), 2026-09-30 02:09:18. Status: **complete (layer probes + head reference)**.

Cells: 90261 of 90261 with an embedding; site4 NK 1690, T 7015; train NK 5640, T 16893.
Look-alike pairs (k = 10): all 1506, no_gdT158 453.

## NK vs T linear probe (site4)

| depth | AUC | balanced acc. | AUC without gdT CD158b+ | pair order rate (all) | pair order rate (no_gdT158) |
|---|---|---|---|---|---|
| in | 0.992 | 0.967 | 0.999 | 0.846 | 0.932 |
| 1 | 0.992 | 0.966 | 0.999 | 0.843 | 0.936 |
| 2 | 0.991 | 0.962 | 0.999 | 0.844 | 0.932 |
| 3 | 0.992 | 0.962 | 0.998 | 0.858 | 0.923 |
| 4 | 0.992 | 0.963 | 0.998 | 0.839 | 0.927 |
| 5 | 0.993 | 0.965 | 0.999 | 0.855 | 0.912 |
| 6 | 0.992 | 0.963 | 0.998 | 0.846 | 0.898 |
| 7 | 0.992 | 0.962 | 0.998 | 0.841 | 0.901 |
| 8 | 0.991 | 0.959 | 0.998 | 0.831 | 0.901 |
| 9 | 0.991 | 0.953 | 0.997 | 0.825 | 0.909 |
| 10 | 0.993 | 0.964 | 0.998 | 0.867 | 0.934 |
| 11 | 0.995 | 0.972 | 0.999 | 0.886 | 0.918 |
| 12 | 0.997 | 0.969 | 0.999 | 0.912 | 0.942 |
| z | 0.997 | 0.969 | 0.999 | 0.920 | 0.938 |
| head outputs (134 predicted proteins) | 0.981 | 0.962 | 0.999 | 0.750 | 0.934 |
| measured proteins (ceiling) | 0.999 | 0.976 | 1.000 | 0.962 | 0.991 |

## Protein probes: site4 Pearson (p95-scaled, clipped)

| depth | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|
| in | 0.695 | 0.670 | 0.518 | 0.865 |
| 1 | 0.688 | 0.665 | 0.533 | 0.866 |
| 2 | 0.680 | 0.664 | 0.530 | 0.857 |
| 3 | 0.675 | 0.662 | 0.521 | 0.856 |
| 4 | 0.676 | 0.657 | 0.513 | 0.857 |
| 5 | 0.680 | 0.666 | 0.521 | 0.863 |
| 6 | 0.676 | 0.664 | 0.527 | 0.860 |
| 7 | 0.674 | 0.663 | 0.522 | 0.859 |
| 8 | 0.672 | 0.662 | 0.516 | 0.855 |
| 9 | 0.659 | 0.657 | 0.507 | 0.846 |
| 10 | 0.687 | 0.659 | 0.537 | 0.858 |
| 11 | 0.682 | 0.663 | 0.528 | 0.870 |
| 12 | 0.717 | 0.685 | 0.582 | 0.896 |
| z | 0.718 | 0.687 | 0.582 | 0.896 |
| head | 0.693 | 0.681 | 0.564 | 0.890 |

## Gap ratio on the look-alike pairs: all (1506 pairs)

Median predicted gap / median measured gap; [95% bootstrap interval]; random NK x T pairs in the last columns.

| depth | CD56 | CD94 | CD335 | CD3 | CD56 random | CD94 random | CD335 random | CD3 random |
|---|---|---|---|---|---|---|---|---|
| in | 0.436 [0.392, 0.525] | 0.293 [0.245, 0.414] | 0.216 [0.195, 0.241] | 0.653 [0.577, 0.729] | 0.637 | 0.753 | 0.551 | 0.506 |
| 1 | 0.416 [0.372, 0.503] | 0.290 [0.239, 0.412] | 0.222 [0.203, 0.248] | 0.641 [0.571, 0.729] | 0.640 | 0.752 | 0.578 | 0.511 |
| 2 | 0.429 [0.386, 0.520] | 0.286 [0.235, 0.404] | 0.216 [0.197, 0.240] | 0.656 [0.588, 0.734] | 0.646 | 0.792 | 0.591 | 0.515 |
| 3 | 0.439 [0.393, 0.528] | 0.302 [0.249, 0.427] | 0.210 [0.191, 0.233] | 0.656 [0.579, 0.732] | 0.646 | 0.784 | 0.572 | 0.517 |
| 4 | 0.421 [0.374, 0.508] | 0.293 [0.244, 0.415] | 0.214 [0.196, 0.238] | 0.668 [0.593, 0.741] | 0.657 | 0.805 | 0.579 | 0.531 |
| 5 | 0.512 [0.447, 0.615] | 0.350 [0.287, 0.506] | 0.245 [0.223, 0.270] | 0.709 [0.632, 0.795] | 0.779 | 0.979 | 0.688 | 0.568 |
| 6 | 0.493 [0.439, 0.597] | 0.356 [0.293, 0.507] | 0.256 [0.234, 0.282] | 0.646 [0.573, 0.720] | 0.780 | 0.977 | 0.709 | 0.563 |
| 7 | 0.459 [0.410, 0.554] | 0.331 [0.273, 0.467] | 0.241 [0.219, 0.270] | 0.649 [0.578, 0.723] | 0.760 | 0.955 | 0.689 | 0.555 |
| 8 | 0.447 [0.398, 0.539] | 0.330 [0.273, 0.471] | 0.225 [0.206, 0.249] | 0.631 [0.561, 0.702] | 0.755 | 0.950 | 0.674 | 0.544 |
| 9 | 0.400 [0.356, 0.488] | 0.294 [0.244, 0.422] | 0.210 [0.190, 0.234] | 0.599 [0.536, 0.677] | 0.708 | 0.917 | 0.632 | 0.509 |
| 10 | 0.462 [0.416, 0.551] | 0.318 [0.262, 0.448] | 0.252 [0.232, 0.276] | 0.671 [0.599, 0.749] | 0.733 | 0.920 | 0.654 | 0.541 |
| 11 | 0.457 [0.405, 0.548] | 0.308 [0.253, 0.439] | 0.259 [0.238, 0.283] | 0.638 [0.572, 0.714] | 0.785 | 0.986 | 0.701 | 0.590 |
| 12 | 0.560 [0.502, 0.662] | 0.337 [0.280, 0.481] | 0.333 [0.304, 0.365] | 0.806 [0.717, 0.896] | 0.955 | 1.152 | 0.930 | 0.729 |
| z | 0.575 [0.520, 0.687] | 0.360 [0.299, 0.511] | 0.340 [0.315, 0.374] | 0.799 [0.718, 0.893] | 0.960 | 1.164 | 0.930 | 0.728 |
| head | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.115 [0.078, 0.147] | 0.154 [0.137, 0.179] | 1.181 | 1.523 | 1.427 | 0.767 |

Measured median gap: CD56 0.189, CD94 0.221, CD335 0.395, CD3 0.137.

## Gap ratio on the look-alike pairs: no_gdT158 (453 pairs)

Median predicted gap / median measured gap; [95% bootstrap interval]; random NK x T pairs in the last columns.

| depth | CD56 | CD94 | CD335 | CD3 | CD56 random | CD94 random | CD335 random | CD3 random |
|---|---|---|---|---|---|---|---|---|
| in | 0.218 [0.182, 0.254] | 0.302 [0.254, 0.386] | 0.293 [0.241, 0.331] | 0.207 [0.188, 0.236] | 0.644 | 0.761 | 0.568 | 0.515 |
| 1 | 0.205 [0.169, 0.246] | 0.280 [0.218, 0.363] | 0.278 [0.236, 0.319] | 0.219 [0.192, 0.243] | 0.649 | 0.763 | 0.596 | 0.521 |
| 2 | 0.205 [0.169, 0.239] | 0.283 [0.232, 0.359] | 0.273 [0.225, 0.314] | 0.209 [0.179, 0.240] | 0.656 | 0.807 | 0.607 | 0.526 |
| 3 | 0.185 [0.153, 0.219] | 0.278 [0.225, 0.354] | 0.241 [0.198, 0.280] | 0.203 [0.179, 0.230] | 0.655 | 0.799 | 0.589 | 0.528 |
| 4 | 0.180 [0.152, 0.213] | 0.269 [0.219, 0.347] | 0.252 [0.198, 0.291] | 0.222 [0.190, 0.246] | 0.668 | 0.824 | 0.597 | 0.543 |
| 5 | 0.198 [0.164, 0.227] | 0.277 [0.228, 0.364] | 0.259 [0.214, 0.289] | 0.224 [0.196, 0.258] | 0.790 | 1.003 | 0.708 | 0.578 |
| 6 | 0.179 [0.154, 0.218] | 0.296 [0.239, 0.372] | 0.255 [0.218, 0.291] | 0.213 [0.181, 0.240] | 0.791 | 0.996 | 0.731 | 0.574 |
| 7 | 0.184 [0.156, 0.215] | 0.292 [0.233, 0.376] | 0.259 [0.215, 0.295] | 0.205 [0.182, 0.236] | 0.772 | 0.973 | 0.708 | 0.566 |
| 8 | 0.172 [0.147, 0.200] | 0.270 [0.215, 0.354] | 0.234 [0.195, 0.274] | 0.197 [0.170, 0.221] | 0.766 | 0.965 | 0.691 | 0.554 |
| 9 | 0.155 [0.128, 0.191] | 0.231 [0.190, 0.305] | 0.217 [0.176, 0.254] | 0.188 [0.160, 0.218] | 0.720 | 0.927 | 0.650 | 0.522 |
| 10 | 0.198 [0.164, 0.237] | 0.287 [0.235, 0.368] | 0.273 [0.229, 0.316] | 0.218 [0.188, 0.244] | 0.745 | 0.932 | 0.667 | 0.553 |
| 11 | 0.213 [0.181, 0.245] | 0.328 [0.270, 0.423] | 0.288 [0.247, 0.328] | 0.221 [0.187, 0.240] | 0.796 | 1.011 | 0.717 | 0.599 |
| 12 | 0.265 [0.228, 0.301] | 0.365 [0.307, 0.468] | 0.390 [0.336, 0.446] | 0.264 [0.224, 0.293] | 0.961 | 1.161 | 0.947 | 0.732 |
| z | 0.279 [0.239, 0.316] | 0.397 [0.327, 0.498] | 0.404 [0.343, 0.460] | 0.276 [0.231, 0.304] | 0.965 | 1.170 | 0.947 | 0.731 |
| head | 0.238 [0.167, 0.312] | 0.311 [0.199, 0.427] | 0.880 [0.695, 1.009] | 0.271 [0.199, 0.364] | 1.180 | 1.528 | 1.438 | 0.768 |

Measured median gap: CD56 0.621, CD94 0.321, CD335 0.413, CD3 0.756.

## Probe minus head (gap ratio, paired bootstrap)

| pair set | depth | CD56 | CD94 | CD335 | CD3 |
|---|---|---|---|---|---|
| all | in | +0.436 [+0.392, +0.525] | +0.293 [+0.245, +0.414] | +0.101 [+0.069, +0.140] | +0.498 [+0.424, +0.574] |
| all | 1 | +0.416 [+0.372, +0.503] | +0.290 [+0.239, +0.412] | +0.107 [+0.075, +0.147] | +0.487 [+0.417, +0.570] |
| all | 2 | +0.429 [+0.386, +0.520] | +0.286 [+0.235, +0.404] | +0.101 [+0.069, +0.141] | +0.502 [+0.433, +0.576] |
| all | 3 | +0.439 [+0.393, +0.528] | +0.302 [+0.249, +0.427] | +0.095 [+0.063, +0.136] | +0.502 [+0.429, +0.576] |
| all | 4 | +0.421 [+0.374, +0.508] | +0.293 [+0.244, +0.415] | +0.099 [+0.069, +0.141] | +0.513 [+0.443, +0.583] |
| all | 5 | +0.512 [+0.447, +0.615] | +0.350 [+0.287, +0.506] | +0.130 [+0.092, +0.174] | +0.554 [+0.483, +0.635] |
| all | 6 | +0.493 [+0.439, +0.597] | +0.356 [+0.293, +0.507] | +0.141 [+0.106, +0.183] | +0.492 [+0.424, +0.561] |
| all | 7 | +0.459 [+0.410, +0.554] | +0.331 [+0.273, +0.467] | +0.126 [+0.093, +0.169] | +0.494 [+0.429, +0.562] |
| all | 8 | +0.447 [+0.398, +0.539] | +0.330 [+0.273, +0.471] | +0.110 [+0.076, +0.150] | +0.477 [+0.407, +0.543] |
| all | 9 | +0.400 [+0.356, +0.488] | +0.294 [+0.244, +0.422] | +0.095 [+0.060, +0.137] | +0.445 [+0.382, +0.517] |
| all | 10 | +0.462 [+0.416, +0.551] | +0.318 [+0.262, +0.448] | +0.137 [+0.102, +0.180] | +0.517 [+0.446, +0.587] |
| all | 11 | +0.457 [+0.405, +0.548] | +0.308 [+0.253, +0.439] | +0.144 [+0.108, +0.187] | +0.484 [+0.419, +0.556] |
| all | 12 | +0.560 [+0.502, +0.662] | +0.337 [+0.280, +0.481] | +0.218 [+0.178, +0.261] | +0.652 [+0.559, +0.738] |
| all | z | +0.575 [+0.520, +0.687] | +0.360 [+0.299, +0.511] | +0.225 [+0.189, +0.269] | +0.645 [+0.563, +0.735] |
| no_gdT158 | in | -0.019 [-0.084, +0.048] | -0.009 [-0.085, +0.112] | -0.587 [-0.713, -0.416] | -0.064 [-0.159, +0.008] |
| no_gdT158 | 1 | -0.032 [-0.096, +0.031] | -0.032 [-0.119, +0.084] | -0.601 [-0.721, -0.429] | -0.052 [-0.144, +0.018] |
| no_gdT158 | 2 | -0.033 [-0.103, +0.031] | -0.028 [-0.106, +0.089] | -0.607 [-0.723, -0.430] | -0.062 [-0.149, +0.005] |
| no_gdT158 | 3 | -0.053 [-0.115, +0.008] | -0.033 [-0.114, +0.084] | -0.639 [-0.760, -0.466] | -0.069 [-0.158, -0.000] |
| no_gdT158 | 4 | -0.058 [-0.119, +0.007] | -0.042 [-0.120, +0.073] | -0.628 [-0.744, -0.459] | -0.050 [-0.141, +0.016] |
| no_gdT158 | 5 | -0.040 [-0.103, +0.023] | -0.034 [-0.108, +0.088] | -0.621 [-0.739, -0.453] | -0.047 [-0.135, +0.025] |
| no_gdT158 | 6 | -0.059 [-0.122, +0.005] | -0.015 [-0.098, +0.107] | -0.625 [-0.733, -0.449] | -0.059 [-0.149, +0.006] |
| no_gdT158 | 7 | -0.054 [-0.118, +0.015] | -0.019 [-0.103, +0.102] | -0.621 [-0.736, -0.451] | -0.066 [-0.153, +0.002] |
| no_gdT158 | 8 | -0.066 [-0.131, -0.001] | -0.041 [-0.123, +0.087] | -0.646 [-0.756, -0.473] | -0.074 [-0.165, -0.006] |
| no_gdT158 | 9 | -0.083 [-0.147, -0.016] | -0.081 [-0.163, +0.043] | -0.663 [-0.784, -0.490] | -0.083 [-0.173, -0.016] |
| no_gdT158 | 10 | -0.040 [-0.107, +0.029] | -0.024 [-0.111, +0.107] | -0.607 [-0.727, -0.432] | -0.053 [-0.147, +0.017] |
| no_gdT158 | 11 | -0.024 [-0.084, +0.035] | +0.017 [-0.060, +0.132] | -0.592 [-0.705, -0.425] | -0.050 [-0.143, +0.014] |
| no_gdT158 | 12 | +0.027 [-0.035, +0.088] | +0.054 [-0.019, +0.179] | -0.489 [-0.597, -0.324] | -0.007 [-0.103, +0.055] |
| no_gdT158 | z | +0.041 [-0.020, +0.102] | +0.086 [+0.002, +0.201] | -0.476 [-0.585, -0.309] | +0.005 [-0.091, +0.068] |

## Pair cosine

| depth | raw: neighbour (all) | raw: random | raw: distance ratio | centred: neighbour (all) | centred: random | centred: distance ratio | own k-NN: NK-T share |
|---|---|---|---|---|---|---|---|
| in | 0.941 | 0.854 | 0.401 | 0.858 | 0.726 | 0.516 | 0.0757 |
| 1 | 0.959 | 0.922 | 0.521 | 0.786 | 0.630 | 0.579 | 0.0824 |
| 2 | 0.972 | 0.948 | 0.534 | 0.712 | 0.500 | 0.576 | 0.0792 |
| 3 | 0.970 | 0.949 | 0.587 | 0.665 | 0.470 | 0.633 | 0.0747 |
| 4 | 0.974 | 0.957 | 0.595 | 0.637 | 0.431 | 0.639 | 0.0697 |
| 5 | 0.979 | 0.961 | 0.543 | 0.630 | 0.403 | 0.620 | 0.0737 |
| 6 | 0.981 | 0.966 | 0.567 | 0.623 | 0.409 | 0.638 | 0.0731 |
| 7 | 0.982 | 0.967 | 0.540 | 0.617 | 0.379 | 0.616 | 0.0753 |
| 8 | 0.980 | 0.964 | 0.537 | 0.612 | 0.364 | 0.610 | 0.0716 |
| 9 | 0.979 | 0.913 | 0.243 | 0.837 | 0.186 | 0.201 | 0.0901 |
| 10 | 0.985 | 0.923 | 0.192 | 0.876 | 0.139 | 0.144 | 0.0728 |
| 11 | 0.986 | 0.914 | 0.164 | 0.878 | 0.240 | 0.161 | 0.0437 |
| 12 | 0.966 | 0.785 | 0.161 | 0.885 | 0.286 | 0.161 | 0.0136 |
| z | 0.966 | 0.785 | 0.161 | 0.883 | 0.274 | 0.161 | 0.0136 |

Ridge alpha at the edge of its grid (the probe may be under- or over-regularised there): in: CD56, CD94, CD3; 1: CD94, CD335, CD3; 2: CD94, CD335, CD3; 3: CD335, CD3; 4: CD94, CD3; 5: CD56, CD3; 6: CD56; 7: CD56, CD3; 8: CD94, CD3; 9: CD3; 10: CD3; 11: CD56, CD3; 12: CD3; z: CD335.

Logistic-regression C at the edge of its grid [0.01, 0.1, 1.0] at depth: in, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11.

Figures: `fig_nkt_probe.svg`, `fig_protein_pearson.svg`, `fig_gap_ratio_all.svg`, `fig_gap_ratio_no_gdT158.svg`, `fig_pair_cosine_raw.svg`, `fig_pair_cosine_centred.svg`, `fig_own_knn_nkt_share.svg`.

## How to read this

**The question.** On site4, NK cells and T cells that TEDDY places next to each other (the NK-T
look-alike pairs of Experiment 6) differ clearly in measured CD56, CD94, CD335 and CD3, but our head's
predictions for the two cells of a pair are much closer than the measurements. Either TEDDY's
representation no longer carries the difference, or it does and our head (MLP + NB decoder trained on z)
throws it away. The probes here read the difference straight off each depth with a linear model and
compare that with the head.

**The numbers.**
- *Gap ratio* (per protein, per pair set): median over pairs of |predicted NK - predicted T| divided by
  the median of |measured NK - measured T|, both on the p95-scaled, clipped scale of the Experiment 6
  tables. 1 = the predictions keep the measured gap; 0 = the two cells of a pair get the same prediction;
  above 1 = the read-out exaggerates the gap (possible where the measured gap is small).
  The interval is a bootstrap over NK cells (each resample keeps all pairs of the drawn NK cells); T cells
  shared between pairs are not clustered, so the interval is on the narrow side.
- *Probe minus head*: probe gap ratio at a depth minus the head's gap ratio, on the same pairs, with a
  paired bootstrap interval (same resamples for both).
- *NK vs T probe*: site4 AUC / balanced accuracy of a logistic regression trained on training-site
  cells; *pair order rate* = share of look-alike pairs in which the NK cell gets the higher P(NK)
  (0.5 = the probe cannot tell the two cells of a pair apart).
- *Cosine*: median cosine of the look-alike pairs and of random NK x T pairs at each depth; the
  distance ratio (1 - cos neighbour) / (1 - cos random) is below 1 when the pairs are closer than chance.

**Reading rule.** Take the gap ratios on the look-alike pairs (both pair sets) and compare the linear
probe with the head:
1. *Loss in the head*: at z (or at an earlier depth) the linear probe's gap ratio is clearly above the
   head's (the probe-minus-head interval lies above 0) and the NK vs T probe orders the look-alike pairs
   well (order rate well above 0.5), while the head's predictions compress the gap. The information is in
   what the head is given; the head does not use it.
2. *Loss in the representation*: the linear probe on z compresses the gap as much as the head does, or
   more (the interval includes 0 or lies below it), and the order rate at z is near 0.5. Then read the
   curve backwards: a depth where the probe still keeps the gap and a later depth where it no longer does
   is where TEDDY (as pooled here) loses it. If no depth keeps it, the difference is not linearly readable
   from any gene-mean layer, including the input embedding, which is only the set of the cell's top genes.
3. *Mixed*: proteins or pair sets can disagree (for example CD3 kept, CD56 lost). Report per protein;
   do not average across them.

**What this does not show.**
- A linear probe is a lower bound on what a depth contains. "Probe compresses too" means the difference
  is not *linearly* available in the gene-mean of that depth; a non-linear read-out or the per-token
  states (not saved) could still hold it. "Probe keeps it" is the stronger statement.
- Ridge shrinks predictions toward the mean, so any probe with imperfect fit compresses gaps. That is why
  the probe is compared with the head on the same pairs, and why the random NK x T pairs are given: a
  look-alike ratio far below the random-pair ratio for the same read-out means the compression is
  specific to look-alikes, not a general loss of the NK-T axis.
- The probes and the head are trained differently (probes: squared error on the p95-scaled, clipped
  scale, one protein at a time; head: NB likelihood on counts for all proteins with a constant size
  factor); both are scored on the same scale and the same pairs.
- The pairs are chosen as neighbours at z, so their cosine at z is high by construction; at earlier
  depths the pairs are not selected and their cosine is only meaningful next to the random-pair cosine.
- Layer L and z are the same states (z is float32 and L2-normalised as the head sees it; layer means
  are stored in float16); the two columns should agree closely, which is a check, not a finding.
- The embedding run used the settings recorded under `inputs.embedding_manifest` (e.g. fp16 autocast);
  the input embedding is recomputed here in float32 from the same tokens.
- One held-out site, one head seed, one pair construction (k = 10). The C and alpha of the probes are
  picked on the val split (a training site), never on site4.

