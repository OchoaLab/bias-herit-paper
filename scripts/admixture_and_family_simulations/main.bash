#!/bin/bash

# admixture and family simulations

m=500000
n=5000
g=1  # g=20 for family simulation
mc=500
t=10
# base name for simulation (structure only)
name=sim-n$n-k3-f0.3-s0.5-g$g
# fuller name for genotypes and phenotypes
name2=$name/m$m

module load R/4.0.0

# create overall population structure params
time Rscript 00-sim-pop.R -n $n -g $g

for rep in {1..50}; do
    # simulate replicate genotypes
    time Rscript 01-sim-gen.R --name $name -r $rep -m $m 
	
    # estimate all kinship matrices to test
    time Rscript 02-kinship-est.R --name $name2 -r $rep -t $t
	
	# simulate replicate phenotypes
	# RC model
	time Rscript 01-sim-phen-tg-rc.R -r $rep --m_causal $mc --herit $h
	
	
	for h in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
	# RC model
	time Rscript 01-sim-phen-rc.R --name $name2 -r $rep --m_causal $mc --herit $h
	# FES model
	# time Rscript 01-sim-phen-fes.R --name $name2 -r $rep --m_causal $mc --herit $h
	done
		
	for h in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
    # calculate heritability estimates
    time Rscript 03-herit-est.R --name $name2 -r $rep -t $t --m_causal $mc --herit $h
	done

done

# gather all small rep tables into a single big one
time Rscript 04-table.R --name $name2 --n_rep 50 --n_h 11 --m_causal $mc

module unload R/4.0.0


