# empirical heritability bias analysis

# admixture case only
time Rscript herit-summary-fig-mira-01-estimate.R -n 1000
# 6m7.481s ideapad
# add minimal family structure (just siblings, so it's fast to generate)
time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 2
# 7m30.000s ideapad
time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 4
# 10m49.391s ideapad
time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 6
# 14m14.933s ideapad
time Rscript herit-summary-fig-mira-01-estimate.R -n 1000 -g 10
# 21m6.778s ideapad
# viiiaX6 has more memory, increased n there only
time Rscript herit-summary-fig-mira-01-estimate.R -n 2000
# 47m16.487s viiiaX6
time Rscript herit-summary-fig-mira-01-estimate.R -n 2000 -g 2
# 39m28.551s viiiaX6
# n=3000 died in viiiaX6, moved to labbyDuke for rest
time Rscript herit-summary-fig-mira-01-estimate.R -n 3000 
# 48m2.400s labbyDuke
time Rscript herit-summary-fig-mira-01-estimate.R -n 3000 -g 2
# 50m7.732s labbyDuke
time Rscript herit-summary-fig-mira-01-estimate.R -n 4000 
# 93m9.443s labbyDuke
time Rscript herit-summary-fig-mira-01-estimate.R -n 4000 -g 2
# 95m43.280s labbyDuke
time Rscript herit-summary-fig-mira-01-estimate.R -n 5000 
# 177m44.662s labbyDuke
time Rscript herit-summary-fig-mira-01-estimate.R -n 5000 -g 2
# 161m4.862s labbyDuke


# plot versions take the same arguments
time Rscript herit-summary-fig-mira-02-plot.R -n 1000 -g 2
# and version of plot actually used in MIRA proposal
time Rscript herit-summary-fig-mira-03-plot-clean.R -n 5000 -g 2


### THEORY ###

# this is for PCA GWAS case

# empirically test the hypothesis that centering and EVD approximation commute
# this script simulates kinship as usual, we use r = K
# RESULTS:
# - commutability is not exact, especially noticeable under family structure, but still appears to be a very good approximation
# - the intercept is in the rowspace in both cases, as expected
# - rowspaces are a mess (centering does alter rowspaces enough that they really don't match; see last line of each run below, calculated via matrix ranks on concatenations of matrices)
# 
# 1) version without family structure (truly low-dimensional, K=3, Fst=0.3)
Rscript test-00-pca-kinship-eqs.R
# RMSD center/dim-r, vs dim-r/center: 0.000501087785158289
# Corr center/dim-r, vs dim-r/center: 0.999993595141344
# Sum of kinship, center/dim-r: 1.3e-14
# Sum of kinship, dim-r/center: -3.7e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank dim-r/center + dim-r: 4
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7
#
# 2) larger K, lower FST, same result
Rscript test-00-pca-kinship-eqs.R -k 10 -f 0.1
# RMSD center/dim-r, vs dim-r/center: 0.0005038297388324
# Corr center/dim-r, vs dim-r/center: 0.999923201852988
# Sum of kinship, center/dim-r: -1.95e-15
# Sum of kinship, dim-r/center: 1.9e-11
# Rank center/dim-r: 10
# Rank dim-r/center: 10
# Rank center/dim-r + dim-r/center: 19
# Rank center/dim-r + dim-r: 20
# Rank dim-r: 10
# Rank dim-r/center: 10
# Rank dim-r/center + dim-r: 11
# Rank intercept + dim-r: 11
# Rank intercept + dim-r/center: 11
# Rank intercept + center/dim-r: 11
# Rank intercept + dim-r/center + dim-r: 11
# Rank intercept + center/dim-r + dim-r: 20
#
# 3) back to default (K=3, Fst=0.3), but add family structure
Rscript test-00-pca-kinship-eqs.R -g 20
# RMSD center/dim-r, vs dim-r/center: 0.00309263509709354
# Corr center/dim-r, vs dim-r/center: 0.99962525703413
# Sum of kinship, center/dim-r: -3.28e-14
# Sum of kinship, dim-r/center: 4.45e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7

