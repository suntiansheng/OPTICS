#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=vary_m_dependent_parallel
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=24:00:00
#SBATCH --output=out/vary_m_dependent_parallel.out
#SBATCH --error=err/vary_m_dependent_parallel.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load gcc/13.3.0
module load openblas/0.3.28
module load r/4.4.1

mkdir -p out err ./results

Rscript -e "install.packages('OPTICS_0.1.1.tar.gz', repos = NULL, type = 'source')"

Rscript vary_m_dependent_simu_parallel.R
