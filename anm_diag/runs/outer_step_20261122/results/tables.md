### Bridge reading: the repaired output, the round against package 2

| target | reference at Z | output in the round (reader of types) | output with the type from the record | beyond the reference | grade |
|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.032 [−0.044, 0.100] support | 0.043 [−0.021, 0.105] unresolved | 0.037 [−0.034, 0.101] unresolved | 0.005 [−0.021, 0.024] no added loss | relative |
| PD-L1, melanoma, co-culture | 0.117 [0.007, 0.193] unresolved | 0.118 [0.020, 0.187] unresolved | 0.108 [0.009, 0.177] unresolved | −0.008 [−0.027, 0.010] no added loss | not tested (gate 3) |
| PD-L1, THP-1 (five events) | 0.132 [0.016, 0.312] unresolved | 0.163 [0.079, 0.290] unresolved | 0.122 [0.000, 0.262] unresolved | −0.010 [−0.064, 0.024] no added loss | relative |
| MHC-I, melanoma, co-culture | 0.321 [0.265, 0.369] reject | 0.391 [0.342, 0.439] reject | 0.362 [0.315, 0.407] reject | 0.040 [0.025, 0.058] no added loss | relative |

### Gate 3 as amended, and the test with two controls

| target | gate 3: output − reference on fitted cells | passes | the output as trained at the same gate | rows − controls-only (gain) | reading | test |
|---|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.03 [−0.09, 0.09] | yes | 0.87 [0.80, 0.91]: fails | 0.080 [0.019, 0.085] | no step | unresolved |
| PD-L1, melanoma, co-culture | −0.09 [−0.15, 0.13] | **no** | 0.88 [0.80, 0.98]: fails | 0.165 [0.150, 0.182] | a step | not testable |
| PD-L1, THP-1 (five events) | −0.12 [−0.16, 0.09] | yes | 0.83 [0.65, 0.94]: fails | 0.372 [0.338, 0.416] | a step | unresolved |
| MHC-I, melanoma, co-culture | −0.08 [−0.09, −0.06] | yes | 0.65 [0.60, 0.70]: fails | 0.276 [0.250, 0.304] | a step | confirmed |
| PD-L1, THP-1 without IRF1 | −0.06 [−0.12, 0.13] | **no** | 0.89 [0.82, 0.95]: fails | 0.364 [0.322, 0.402] | a step | not testable |

### What each arm removes (gain: fall of \|*u*\| of the output)

| target | as trained | controls only | rows | every cell | rows of the other group | the other group's repair |
|---|---|---|---|---|---|---|
| PD-L1, melanoma, IFN-γ | 0.91 | 0.79 (leaves 0.12) | 0.87 (leaves 0.04) | 0.90 (leaves 0.00) | 0.89 (leaves 0.02) | 0.63 (leaves 0.28) |
| PD-L1, melanoma, co-culture | 1.00 | 0.72 (leaves 0.27) | 0.89 (leaves 0.11) | 0.89 (leaves 0.11) | 0.84 (leaves 0.16) | 0.65 (leaves 0.35) |
| PD-L1, THP-1 (five events) | 0.97 | 0.47 (leaves 0.49) | 0.84 (leaves 0.12) | 0.86 (leaves 0.10) | 0.73 (leaves 0.23) | 0.21 (leaves 0.76) |
| MHC-I, melanoma, co-culture | 0.97 | 0.33 (leaves 0.64) | 0.61 (leaves 0.36) | 0.63 (leaves 0.34) |  |  |

### The THP-1 pool, three ways

| read as | reference | repaired output | beyond the reference |
|---|---|---|---|
| PD-L1, THP-1 (five events) | 0.13 [0.02, 0.31] unresolved | 0.12 [0.00, 0.26] unresolved | −0.01 [−0.06, 0.02] no added loss |
| PD-L1, THP-1 without IRF1 | 0.08 [−0.03, 0.15] unresolved | 0.09 [−0.04, 0.16] unresolved | 0.01 [−0.02, 0.04] no added loss |
| IRF1 with its own row | 0.40 [0.27, 0.50] reject (a reader of Z with its own cells) | 0.68 [0.57, 0.79] reject | 0.28 [0.20, 0.39] loss |
| IRF1 with the pooled rows | the same | 0.54 [0.43, 0.65] reject | 0.15 [0.08, 0.22] unresolved |
| IRF1, every cell of the fitted half | the same | 0.42 [0.30, 0.53] reject | 0.02 [−0.04, 0.09] no added loss |

IRF1 before the head: the RNA view (512 components, ridge, own cells) leaves 0.18 [0.04, 0.30] unresolved; Z leaves 0.40 [0.27, 0.50]; from the RNA view to Z 0.21 [0.12, 0.32]: loss.

### Readers: PD-L1 in co-culture (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| IFNGR1 | 0.10 unr | 0.18 unr | 0.19 rej | 0.09 unr | 0.13 unr | 1.00 | 0.01 [−0.04, 0.05] |
| IFNGR2 | 0.11 unr | 0.17 unr | 0.21 rej | 0.14 unr | 0.16 unr | 1.00 | −0.04 [−0.08, 0.01] |
| JAK1 | −0.06 unr | −0.06 unr | 0.04 unr | −0.07 unr | −0.07 unr | 1.00 | 0.02 [−0.11, 0.15] |
| JAK2 | 0.14 unr | 0.18 rej | 0.23 rej | 0.15 unr | 0.19 rej | 1.00 | −0.01 [−0.05, 0.03] |
| STAT1 | 0.18 rej | 0.21 rej | 0.23 rej | 0.16 unr | 0.19 rej | 1.00 | 0.01 [−0.03, 0.05] |
| pool | 0.10 unr | 0.14 unr | 0.19 rej | 0.11 unr | 0.13 unr | 1.00 | −0.00 [−0.04, 0.03] |

Positive control: of 60 readings (two scopes, five readers, every event and the pool) 0 supported, 40 unresolved, 20 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### Readers: MHC-I in co-culture (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| IFNGR1 | 0.25 rej | 0.29 rej | 0.27 rej | 0.13 unr | 0.27 rej | 1.00 | 0.12 [0.07, 0.16] |
| IFNGR2 | 0.32 rej | 0.37 rej | 0.33 rej | 0.22 rej | 0.34 rej | 1.00 | 0.10 [0.07, 0.14] |
| JAK1 | 0.34 rej | 0.39 rej | 0.32 rej | 0.20 rej | 0.32 rej | 1.00 | 0.13 [0.10, 0.17] |
| JAK2 | 0.29 rej | 0.36 rej | 0.30 rej | 0.18 rej | 0.31 rej | 1.00 | 0.11 [0.07, 0.15] |
| STAT1 | 0.35 rej | 0.39 rej | 0.33 rej | 0.22 rej | 0.35 rej | 1.00 | 0.12 [0.09, 0.16] |
| pool | 0.31 rej | 0.37 rej | 0.31 rej | 0.19 rej | 0.32 rej | 1.00 | 0.12 [0.09, 0.15] |

Positive control: of 60 readings (two scopes, five readers, every event and the pool) 0 supported, 2 unresolved, 58 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### Readers: MHC-I under IFN-γ: B2M and TAPBP (own cells)

| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |
|---|---|---|---|---|---|---|---|
| B2M | 0.91 rej | 0.97 rej | 1.00 rej | 1.00 rej | 0.99 rej | 1.00 | −0.07 [−0.10, −0.04] dominated |
| TAPBP | 1.00 rej | 1.00 rej | 1.02 rej | 1.03 rej | 1.04 rej | 1.00 | −0.02 [−0.07, 0.02] |
| pool | 0.94 rej | 0.98 rej | 1.00 rej | 1.01 rej | 1.00 rej | 1.00 | −0.06 [−0.09, −0.01] dominated |

Positive control: of 30 readings (two scopes, five readers, every event and the pool) 0 supported, 0 unresolved, 30 rejected: does not hold. Permuted control: holds. The input view cannot be read.

### The selection third

| readout | knockouts that pass A1 on a third | admitted | the next one |
|---|---|---|---|
| PDL1 | CD274, IFNGR1, STAT1, IFNGR2, JAK2, JAK1 | IFNGR1, STAT1, IFNGR2, JAK2, JAK1 | NSG1 (1.5) |
| MHC1 | B2M, TAPBP | B2M, TAPBP | HLA-B (2.3) |
| HLAE | none | none | STAT1 (1.6) |

| events of the round, on a third | |
|---|---|
| PDL1 mel_ifn | IFNGR1 5.6, IFNGR2 4.8, JAK1 3.1, JAK2 4.6, STAT1 5.0 |
| PDL1 mel_co | IFNGR1 5.9, IFNGR2 5.5, JAK1 5.0, JAK2 5.6, STAT1 6.9 |
| MHC1 mel_co | IFNGR1 8.5, IFNGR2 7.6, JAK1 9.4, JAK2 8.2, STAT1 10.9 |
| PDL1 thp1 | IFNGR1 7.6, IFNGR2 8.6, IRF1 2.6, JAK2 7.6, STAT1 4.3 |

### The round under package 2: MHC-I in melanoma under IFN-γ

| event | cells | change (control SD) | head as trained | Z, other cells of the type | Z, own cells | A3 | rows | controls only | every cell | test |
|---|---|---|---|---|---|---|---|---|---|---|
| B2M | 168 | −2.88 | 1.00 [0.99, 1.00] | 1.00 [0.97, 1.02] | 0.91 [0.88, 0.94] | fails (−0.00) | 0.99 [0.97, 1.02] | 1.01 [0.99, 1.03] | 0.92 [0.89, 0.95] | not testable |
| TAPBP | 153 | −1.60 | 1.00 [0.99, 1.02] | 1.01 [0.97, 1.04] | 1.00 [0.96, 1.04] | fails (−0.00) | 1.01 [0.97, 1.06] | 1.01 [0.97, 1.06] | 1.01 [0.97, 1.05] | not testable |

### The census and the final units

| stratum | units (at least 100 cells) | read by the system | unread | final | control cells | with states | unread |
|---|---|---|---|---|---|---|---|
| frangieh2021, Co-culture | 197 | 5 | 192 | 64 | 5897 | 4000 | 1897 |
| frangieh2021, Control | 188 | 0 | 188 | 63 | 3978 | 0 | 3978 |
| frangieh2021, IFNγ | 206 | 206 | 0 | 0 | 6444 | 4000 | 2444 |
| papalexi2021, all | 24 | 5 | 19 | 6 | 2386 | 2386 | 0 |

| stratum | batch 1: chosen for an open item | calibration | screened as clear, but final |
|---|---|---|---|
| frangieh2021, Co-culture | CD44, TAPBP | CDH19 | B2M, NPC1 |
| frangieh2021, Control | B2M | KLF4 | STAT3, TAPBP |
| papalexi2021, all | BRD4, CMTM6 | POU2F2 | CD86 |

### The forecasts

| hypothesis | intervals inside | outcomes: the class given the highest probability happened |
|---|---|---|
| E1 | 0 of 0 | 5 of 5 |
| E2 | 1 of 5 | 6 of 6 |
| E3 | 3 of 3 | 0 of 0 |
| P1 | 2 of 4 | 7 of 8 |
| P2 | 4 of 4 | 6 of 8 |
| P4 | 1 of 1 | 4 of 4 |
| S1 | 2 of 3 | 0 of 2 |
| S2 | 3 of 4 | 4 of 6 |
| V1 | 1 of 2 | 3 of 3 |
| V2 | 2 of 3 | 2 of 3 |
| all | 19 of 29 | 37 of 45 |

| missed | forecast | observed |
|---|---|---|
| P1: gate 3 passes, PDL1 mel_co | yes 0.95, no 0.05 | no |
| P1: gate 3 against the reference, PDL1 thp1 | [-0.11, -0.03] | −0.117 |
| P1: gate 3 against the reference, MHC1 mel_co | [-0.07, -0.01] | −0.079 |
| P2: outcome of the test with two controls, PDL1 mel_co | unresolved 0.65, confirmed 0.35 | not testable |
| P2: outcome of the test with two controls, MHC1 mel_co | unresolved 0.7, confirmed 0.3 | confirmed |
| S1: PDL1 mel_ifn verdict | support 0.6, unresolved 0.4 | unresolved |
| S1: PDL1 thp1 (five events pooled), repaired output | [0.14, 0.17] | 0.122 |
| S1: PDL1 mel_ifn closed (absolute) on development cells | yes 0.55, no 0.45 | no |
| S2: type B (IRF1, its own row) repaired output | [0.35, 0.55] | 0.678 |
| S2: type B adds no loss beyond a reader of Z with its own cells | yes 0.6, no 0.4 | no |
| S2: type B reading upstream | yes 0.45, no 0.55 | yes |
| V1: PDL1 mel_co, 512 components, boosting, own cells | [0.15, 0.3] | 0.107 |
| V2: MHC1 mel_co, 512 components, boosting, own cells | [0.22, 0.42] | 0.194 |
| V2: boosting reads at least 0.10 better than ridge at a stage | yes 0.2, no 0.8 | yes |
| E2: head as trained, TAPBP | [0.85, 1.0] | 1.002 |
| E2: Z, own cells, B2M | [0.3, 0.9] | 0.913 |
| E2: Z, own cells, TAPBP | [0.4, 0.95] | 0.997 |
| E2: input view, own cells, B2M (best reader) | [0.2, 0.85] | 0.986 |
