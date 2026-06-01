# m-dependent

m-dependent noise changepoint simulations.

## Entry Points

- `m_dependent_simu_parallel.R` (primary)
- `m_dependent_simu_parallel_d5.R` (d=5 variant)
- `job.sh`

## Run

```bash
cd m-dependent
sbatch job.sh
# or
Rscript m_dependent_simu_parallel.R
```

Run d=5 variant locally if needed:

```bash
Rscript m_dependent_simu_parallel_d5.R
```

## Outputs

- `results/m_dependent.csv`
- `results/m_dependent_d5.csv` (when d=5 script is run)
