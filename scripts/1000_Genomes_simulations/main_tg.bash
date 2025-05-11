#!/bin/bash

# 1000 Genomes simulations

mc=500
h=0.8
t=10

module load R/4.0.0

# preprocess genotypes for different MAF filters 

# estimate kinship matrices using GCTA and popkin R package 

for rep in {1..50}; do

	# simulate replicate phenotypes
	# RC model
	time Rscript 01-sim-phen-tg-rc.R -r $rep --m_causal $mc --herit $h
	# FES model
	#time Rscript 01-sim-phen-tg-fes.R -r $rep --m_causal $mc --herit $h

	# calculate heritability estimates
	time Rscript 03-herit-est-tg.R -r $rep -t $t --m_causal $mc --herit $h
	
done

# gather all small rep tables into a single big one

module unload R/4.0.0


