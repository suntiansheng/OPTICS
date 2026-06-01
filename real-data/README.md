# real-data

Real-data analysis using the `ACGH` dataset from the `ecp` package.

## Entry Point

- `real_data.R`

## Run

```bash
cd real-data
Rscript real_data.R
```

## Inputs

- `data(ACGH, package = "ecp")` loaded in-script.

## Outputs

Figures in `figure/`:
- `moment_upper_bound.png`
- `moment_upper_bound.pdf`
- `combined_boxplot.png`
- `combined_boxplot.pdf`
- `SNP_line_changepoints1_10.png`
- `SNP_line_changepoints1_10.pdf`
- `joint.png`
- `joint.pdf`

Serialized objects in `result/`:
- `cpt_ls.RData`
- `A_ls.RData`
- `min_ls.RData`
- `A_joint_OPTICS.RData`
- `cp_joint_WBS.RData`
