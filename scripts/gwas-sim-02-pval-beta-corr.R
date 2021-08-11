library(optparse) # for terminal options
library(readr)    # to read tables
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF

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
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model")
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
fes <- opt$fes

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
    if ( fes ) '-fes' else '-rc'
)

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )
n_kinship <- nrow( kinship_methods )

# let's plot data in the order of the `kinship_methods` table
# also, plot PCA first, LMM second
method_codes <- c(
    paste0( 'pca_', kinship_methods$code ),
    paste0( 'lmm_', kinship_methods$code )
)
# same but human-readable
method_nice <- c(
    paste0( 'PCA, ', kinship_methods$nice ),
    paste0( 'LMM, ', kinship_methods$nice )
)

# load pre-existing data
setwd( '../data/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/true-vs-biased-kinship-gwas/results' )
setwd( dir_out )

# load tibbles
pvals <- read_tsv( 'pvals.txt', col_types = cols( ) )
betas <- read_tsv( 'betas.txt', col_types = cols( ) )

plot_cor <- function( data, name ) {
    # reorder columns of data to be a manually-selected order
    data <- data[ method_codes ]
    # replace original codes with nice names
    colnames( data ) <- method_nice
    # compute correlation matrix
    cor_data <- cor( data )
    
    # get nice max width for a journal
    dim <- fig_width()
    fig_start(
        name,
        width = dim,
        height = dim * 0.84
    )
    plot_popkin(
        cor_data,
        names = TRUE,
        names_cex = 0.7,
        mar = 8,
        ylab = 'Association Model, Kinship Estimate',
        ylab_adj = 0.75,
        leg_title = expression(bold(paste("Pearson Correlation ", (rho)))),
        leg_width = 0.15
    )
    fig_end()

    # return in same order, etc
    return( cor_data )
}

cor_pvals <- plot_cor( pvals, 'pvals_cor' )
cor_betas <- plot_cor( betas, 'betas_cor' )

# pick out some values of particular importance
# p-values only
range_subset <- function( indexes ) {
    # take subset
    cor_pvals <- cor_pvals[ indexes, indexes ]
    # be clear about what is being included:
    message( 'Subset: ', toString( rownames( cor_pvals ) ) )
    message( 'Range: ', toString( range( cor_pvals ) ) )
}
# these should be all PCA methods except for True and Popkin
range_subset( c( 2:4, 6:9 ) )
# include Popkin too
range_subset( 2:9 )
# include True too
range_subset( 1:9 )
# now LMM subsets
# limits only except GCTA
range_subset( 10:12 )
# all limits only
range_subset( 10:13 )
# estimates ROM: popkin-wg-std
range_subset( 14:16 )
# estimates MOR: std-gcta
range_subset( 17:18 )
# all estimators only
range_subset( 14:18 )
# all LMM
range_subset( 10:18 )
