library(optparse) # for terminal options
library(readr)    # to read tables
library(gplots)   # to plot
#library(popkin)   # to plot
library(ochoalabtools) # for nice PDF
library(RColorBrewer)

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 10000, 
                help = "number of loci", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 3, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.3, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int"),
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
const_herit_loci <- opt$const_herit_loci

# output path for BED files and all results files
dir_out <- paste0(
    'sim-admix',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-g', generations,
    if ( const_herit_loci ) '-inv' else '-rand'
)

# load pre-existing data
setwd( '../data/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/true-vs-biased-kinship-gwas/results' )
setwd( dir_out )

# load tibbles
pvals <- read_tsv( 'pvals.txt' )
betas <- read_tsv( 'betas.txt' )

plot_cor <- function( data, name ) {
    # get nice max width for a journal
    dim <- fig_width()
    fig_start(
        name,
        width = dim,
        height = dim,
        mar_l = 5,
        mar_b = 5
    )
    heatmap.2(
        cor( data ),
        symm = TRUE,
        col = brewer.pal(9, "Reds"),
        #    breaks = seq(-1,1,0.02),
        key.title = NA,
        key.xlab = 'Correlation',
        margins = c(8, 8),
        offsetRow = 0,
        offsetCol = 0,
        trace = "none",
        density.info = "none"
    )
    fig_end()
}

plot_cor( pvals, 'pvals_cor' )
plot_cor( betas, 'betas_cor' )
