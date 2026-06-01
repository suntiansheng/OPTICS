# mean

Mean changepoint simulations under normal and heavy-tail settings.

## Entry Points

- `mean_simu_parallel.R`
- `job.sh`

## Run

```bash
cd mean
sbatch job.sh
# or
Rscript mean_simu_parallel.R
```

## Outputs

- `results/multi_d_norm.csv`
- `results/multi_d_t.csv`
- Optional commented paths in script for 1D outputs:
  - `results/one_d_normal.csv`
  - `results/one_d_t.csv`
