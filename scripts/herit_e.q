#!/bin/bash
#SBATCH -p biostat --account=biostat
#SBATCH --job-name=herit_e%a
#SBATCH --output=herit_e%a.out
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --mail-user=zh105@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 16G was not enough locally

m=500000
n=5000
g=1
mc=500
h=0.8
t=10
# base name for simulation (structure only)
name=sim-n$n-k3-f0.3-s0.5-g$g
# fuller name for genotypes and phenotypes
name2=$name/m$m

rep=$SLURM_ARRAY_TASK_ID

module load R/4.0.0

	for h in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do


    # calculate heritability estimates!
    time Rscript 03-herit-est.R --name $name2 -r $rep -t $t --m_causal $mc --herit $h
    # 34m56.756s/54m31.964s viiiaR5
	done


module unload R/4.0.0
