# multiple

Multiple-splitting checks and baseline-vs-multiple comparisons.

## Entry Points

- `multiple_simu.R`
- `multiple_simu_parallel.R`
- `check_ms.R`

## Run

```bash
cd multiple
Rscript multiple_simu.R
Rscript check_ms.R
```

## Outputs

- `results/mean_multiple_compare_summary.csv`
- `results/check_ms_detail.csv`
- `results/check_ms_summary.csv`
- Optional parallel outputs from `multiple_simu_parallel.R`:
  - `results/one_d_normal_multiple.csv`
  - `results/one_d_t10_multiple.csv`
