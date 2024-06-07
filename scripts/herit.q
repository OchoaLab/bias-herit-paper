#!/bin/bash
#SBATCH -p biostat --account=biostat
#SBATCH --job-name=herit%a
#SBATCH --output=herit%a.out
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

# control threads on DCC (default is all CPUs in node!)
    # simulate replicate genotypes and phenotypes
    # NOTE: old code looped inside but ran out of memory (OOM) consistently, but only after one iteration.  This solution of looping outside of R successfully fixed OOM problems!
    time Rscript 01-sim-gen.R --name $name -r $rep -m $m 
    # 0m54.031s/4m41.821s viiiaR5
	

    # estimate all kinship matrices to test
    time Rscript 02-kinship-est.R --name $name2 -r $rep -t $t
    # 2m32.873s/7m47.149s viiiaR5
	
	for h in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
	time Rscript 01-sim-phen.R --name $name2 -r $rep --m_causal $mc --herit $h
	done



module unload R/4.0.0
