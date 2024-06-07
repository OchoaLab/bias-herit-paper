

The main.bash is for reference. I actually didn't run that directly but ran jobs seperately. 
Here is my workflow.

1.Run this in the terminal (on interactive session, cpu=40, Memory=208)

```bash 

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



```

2.Edit herit.q file and run this:

herit.q: simulate genotypes and phenotypes and estimate all kinship matrices

``` bash

cd /work/zh105/Alex_projects/project2/scripts/


sbatch -a 1-50 herit.q


```
3.Calculate heritability estimates using herit_e.q

```bash
sbatch -a 1-50 herit_e.q


```


4.Summarize all results

```bash

cd /work/zh105/Alex_projects/project2/scripts/

m=500000
n=5000
g=1
mc=500
h=0.8
t=5
# base name for simulation (structure only)
name=sim-n$n-k3-f0.3-s0.5-g$g
# fuller name for genotypes and phenotypes
name2=$name/m$m

module load R/4.0.0


time Rscript 04-table.R --name $name2 --n_rep 50 --n_h 11 --m_causal $mc


```

5.Make plots using the big table:
05-boxplots.R

All the raw code, data and figures are within this folder.


