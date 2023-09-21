# empirical heritability bias analysis

n=5000
g=2
# base name for simulation (structure only)
name=sim-n$n-k3-f0.3-s0.5-g$g
# fuller name for genotypes and phenotypes
name2=$name/m100000-mc100-h0.8

# create overall population structure params
time Rscript 00-sim-pop.R -n $n -g $g
# 0m15.793s viiiaR5

for rep in {1..50}; do
    # simulate replicate genotypes and phenotypes
    # NOTE: old code looped inside but ran out of memory (OOM) consistently, but only after one iteration.  This solution of looping outside of R successfully fixed OOM problems!
    time Rscript 01-sim-gen-phen.R --name $name -r $rep
    # 0m54.031s/4m41.821s viiiaR5

    # estimate all kinship matrices to test
    time Rscript 02-kinship-est.R --name $name2 -r $rep
    # 2m32.873s/7m47.149s viiiaR5

    # calculate heritability estimates!
    time Rscript 03-herit-est.R --name $name2 -r $rep
    # 34m56.756s/54m31.964s viiiaR5
done

# gather all small rep tables into a single big one
time Rscript 04-table.R --name $name2 --n_rep 50
# 0m0.874s viiiaR5

# boxplot of heritability estimates!
# creates smaller "grant" version as well (no separate flag in this case)
time Rscript 05-herit-boxplots.R --name $name2



### OLDER

# left this as a sample DCC submission script, but it's not currently used
sbatch herit.q

# separate analyses

# confirm that kinship estimation noise leads to bias, is reduced to zero as number of loci increase (absent all other sources of bias)
time Rscript test-kinship-noise.R
# 9m9.047s/62m39.499s viiiaR5
time Rscript test-kinship-noise.R -g 2
# not yet run
