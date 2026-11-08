# Synthetic data

`daily_synthetic.csv` is simulated by `R/02_make_synthetic_data.R` (seed 2026).
It has the same columns as the real processed file, but contains no real observations.

- Predictors are drawn from a multivariate normal using summary statistics of the real data.
- Biocharge is built from a KNOWN, invented formula, so models can be checked on whether they recover it.
- The formula is NOT Zepp's. It exists only to test the analysis pipeline.
- 1 to 14 May are missing, as in the real data.
- Claude (Anthropic) assisted in synthesising the data, just to speed up time.
