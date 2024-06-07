# empirical heritability bias analysis


cd /work/zh105/Alex_projects/project2/scripts/


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

# create overall population structure params
time Rscript 00-sim-pop.R -n $n -g $g



for rep in {1..50}; do
    # simulate replicate genotypes and phenotypes
    # NOTE: old code looped inside but ran out of memory (OOM) consistently, but only after one iteration.  This solution of looping outside of R successfully fixed OOM problems!
    time Rscript 01-sim-gen.R --name $name -r $rep -m $m 
    # 0m54.031s/4m41.821s viiiaR5
	

    # estimate all kinship matrices to test
    time Rscript 02-kinship-est.R --name $name2 -r $rep -t $t
    # 2m32.873s/7m47.149s viiiaR5
	
	for h in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
	time Rscript 01-sim-phen.R --name $name2 -r $rep --m_causal $mc --herit $h
	

    # calculate heritability estimates!
    time Rscript 03-herit-est.R --name $name2 -r $rep -t $t --m_causal $mc --herit $h
    # 34m56.756s/54m31.964s viiiaR5
	done
done

# gather all small rep tables into a single big one
time Rscript 04-table.R --name $name2 --n_rep 3 --n_h 11 --m_causal $mc
# 0m0.874s viiiaR5

# boxplot of heritability estimates!
# creates smaller "grant" version as well (no separate flag in this case)
time Rscript 05-boxplots.R --name $name2


