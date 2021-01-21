#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu=2G # 16 GB RAM
module load R/4.0.0
Rscript heritability_MVN_genetic_bias_1221_const_herit_loci.R -t 8 -g 2 -r 10 -n 5000


