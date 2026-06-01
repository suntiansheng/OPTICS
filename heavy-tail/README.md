# heavy-tail

Heavy-tail mean changepoint simulation (t-distributed errors).

## Entry Points

- `heavy_simu_parallel.R`
- `job.sh`

## Run

```bash
cd heavy-tail
sbatch job.sh
# or
Rscript heavy_simu_parallel.R
```

## Output

- `results/heavy_one_d.csv`
- Repository currently also contains curated export `heavy_one_d.csv`.
