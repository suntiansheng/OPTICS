# covariance

Covariance changepoint simulation with Toeplitz alternatives.

## Entry Points

- `covariance_simu_parallel.R`
- `job.sh`

## Run

```bash
cd covariance
sbatch job.sh
# or
Rscript covariance_simu_parallel.R
```

## Output

- `results/covariance.csv`
- Repository currently also contains curated export `covariance.csv`.
