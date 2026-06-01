# linear

Linear signal changepoint simulations for normal and t settings.

## Entry Points

- `linear_simu_parallel.R`
- `job.sh`

## Run

```bash
cd linear
sbatch job.sh
# or
Rscript linear_simu_parallel.R
```

## Outputs

- `results/linear_normal.csv`
- `results/linear_t.csv`
