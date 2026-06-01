# vary-m-dependent

Simulation varying the dependence order parameter `M`.

## Entry Points

- `vary_m_dependent_simu_parallel.R`
- `job.sh`

## Run

```bash
cd vary-m-dependent
sbatch job.sh
# or
Rscript vary_m_dependent_simu_parallel.R
```

## Output

- `results/vary_m_dependent.csv` (script output path)
- Repository currently also contains curated export `vary_m_dependent.csv`.
