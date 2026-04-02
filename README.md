# OPTICS Simulation and Real-Data Experiments

This repository contains simulation studies and real-data analyses for changepoint model selection using `OPTICS` and related methods.

## Repository Layout

- `mean/`: mean changepoint simulations (1D and multi-D variants)
- `var/`: variance changepoint simulation
- `linear/`: linear signal simulation
- `network/`: network-structured simulation
- `m-dependent/`: m-dependent noise simulations (1D and d=5)
- `vary-m-dependent/`: varying m-level simulation
- `covariance/`: covariance changepoint simulation
- `heavy-tail/`: heavy-tail mean changepoint simulation
- `multiple/`: multiple-splitting checks and comparisons
- `MS-OPTICS/`: summarized multiple-splitting benchmark pipeline
- `real-data/`: analysis and figures for `ecp::ACGH` data
- `OPTICS_0.1.1.tar.gz`: local package tarball used by job scripts
- `CODE_TIDY_REQUIREMENTS.md`: code hygiene and reproducibility standards

## Environment

- R >= 4.4 recommended
- Package tarball: `OPTICS_0.1.1.tar.gz`
- Additional R packages used in scripts include `parallel`, `FDRSeg`, `MASS`, `changepoints`, `changepoint`, `ecp`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, and `pracma`.

## Running Experiments

For HPC execution, each simulation folder with a `job.sh` can be run with:

```bash
cd <folder>
sbatch job.sh
```

For local execution:

```bash
cd <folder>
Rscript <script_name>.R
```

## Main Output Locations

- Simulation tables are written under each folder's `results/` directory (or folder root for some curated CSV exports).
- Real-data outputs are written to:
  - `real-data/figure/`
  - `real-data/result/`

## Reproducibility Notes

- Scripts use fixed seeds and explicit model grids.
- Parallel scripts initialize worker RNG streams where applicable.
- Keep output file names and column schemas unchanged across reruns.
- Follow `CODE_TIDY_REQUIREMENTS.md` before final manuscript submission.
