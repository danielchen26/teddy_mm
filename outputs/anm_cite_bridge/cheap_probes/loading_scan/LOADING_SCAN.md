# Loading scan (Mode B)

**Not Mode A.**

## Re-encode feasibility (ctx 256/512/1024)

- feasible: **False**
- blocker: `torch_unavailable: ModuleNotFoundError: No module named 'torch'`
- note: `ckpt_weights_bytes=284887424`
- note: `processed_cite_arrays_present`
- available on disk: `[{'ctx': 1024, 'path': '/Users/tianchichen/Documents/GitHub/teddy_mm/data/processed/cite/z_rna.npy', 'note': 'scripts/03_embed_rna.py default'}]`

## Proxy: leading-dim subsets of ctx1024 z_512

| variant | dim | near | must | frac must |
|---|---:|---:|---:|---:|
| leading_dim_128 | 128 | 116434 | 50003 | 0.4294535960286514 |
| leading_dim_256 | 256 | 155262 | 68780 | 0.44299313418608544 |
| leading_dim_512 | 512 | 177700 | 79048 | 0.444839617332583 |

must_frac_spread: **0.015386021303931574** | rel small vs 512: **-0.03458779457682219**

Do not claim Pearson > ~0.61. Not clinical.
