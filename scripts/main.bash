# empirical heritability bias analysis

# not rerun since last simtrait/simfam updates
# # admixture case only
# time Rscript herit-summary-fig-mira-01-estimate.R -n 1000
# # 6m7.481s ideapad
# # add minimal family structure (just siblings, so it's fast to generate)
# time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 2
# # 7m30.000s ideapad
# time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 4
# # 10m49.391s ideapad
# time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 6
# # 14m14.933s ideapad
# time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 10
# # 21m6.778s ideapad
# # viiiaX6 has more memory, increased n there only
# time Rscript herit-summary-fig-mira-01-estimate.R -n 2000
# # 47m16.487s viiiaX6
# time Rscript herit-summary-fig-mira-01-estimate.R -n 2000 -g 2
# # 39m28.551s viiiaX6
# # n=3000 died in viiiaX6, moved to labbyDuke for rest
# time Rscript herit-summary-fig-mira-01-estimate.R -n 3000 
# # 48m2.400s labbyDuke
# time Rscript herit-summary-fig-mira-01-estimate.R -n 3000 -g 2
# # 50m7.732s labbyDuke
# time Rscript herit-summary-fig-mira-01-estimate.R -n 4000 
# # 93m9.443s labbyDuke
# time Rscript herit-summary-fig-mira-01-estimate.R -n 4000 -g 2
# # 95m43.280s labbyDuke
# time Rscript herit-summary-fig-mira-01-estimate.R -n 5000 
# # 177m44.662s labbyDuke
# time Rscript herit-summary-fig-mira-01-estimate.R -n 5000 -g 2
# # 161m4.862s labbyDuke
sbatch herit.q # same as last above but on DCC
# 1042m54.307s/1156m26.282s DCC

# plot versions take the same arguments
# time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 1
# time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 2
# time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 4
# time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 6
# time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 10
# time Rscript herit-summary-fig-mira-02-plot.R -n 2000 -g 1
# time Rscript herit-summary-fig-mira-02-plot.R -n 2000 -g 2
# time Rscript herit-summary-fig-mira-02-plot.R -n 3000 -g 1
# time Rscript herit-summary-fig-mira-02-plot.R -n 3000 -g 2
# time Rscript herit-summary-fig-mira-02-plot.R -n 4000 -g 1
# time Rscript herit-summary-fig-mira-02-plot.R -n 4000 -g 2
# time Rscript herit-summary-fig-mira-02-plot.R -n 5000 -g 1
time Rscript herit-summary-fig-mira-02-plot.R -n 5000 -g 2
# and version of plot actually used in MIRA proposal
time Rscript herit-summary-fig-mira-02-plot.R -n 5000 -g 2 --grant

# separate analyses

# confirm that kinship estimation noise leads to bias, is reduced to zero as number of loci increase (absent all other sources of bias)
time Rscript test-kinship-noise.R
# 9m9.047s/62m39.499s viiiaR5
time Rscript test-kinship-noise.R -g 2
# not yet run
