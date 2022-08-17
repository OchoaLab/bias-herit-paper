#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=herit
#SBATCH --output=herit.out
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=6
#SBATCH --mail-user=alejandro.ochoa@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 16G was not enough locally

module load R/4.0.0

# control threads on DCC (default is all CPUs in node!)
time Rscript herit-summary-fig-mira-01-estimate.R -n 5000 -g 2 -t 6

module unload R/4.0.0
