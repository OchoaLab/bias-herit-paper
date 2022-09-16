library(optparse) # for terminal options
library(readr)    # to read tables
library(ochoalabtools) # for nice PDF
library(dplyr)    # for bind_rows
library(popkin)
library(genio) # read_grm, to get true mean kinship

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--name", type = "character", default = NA, 
                help = "Base name for genotype and phenotype simulation", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$name
n_rep <- opt$n_rep
if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )
# reorder for desired plotting order
kinship_methods <- kinship_methods[ order( match( kinship_methods$short, c('True Kinship', 'Popkin', 'Standard', 'Weir-Goudet')) ), ]


# go where the data is
setwd( '../data/' )
setwd( dir_out )

# recover true heritability from dir name
herit <- as.numeric( sub( '.*-h', '', dir_out ) )

# also calculate mean kinship for each original mat, to predict bias?
# just focus on the real one?
mean_kinship_true <- mean( read_grm( '../kinship/true' )$kinship / 2 )
# actually predict biased herit value from formula
herit_biased <- herit * ( 1 - mean_kinship_true ) / ( 1 - ( mean_kinship_true * herit ) )

# tibble to grow
data <- NULL

# load pre-calculated herit estimates
for ( rep in 1 : n_rep ) {
    file_rep <- paste0( 'rep-', rep, '/herit.txt.gz' )
    data_rep <- read_tsv( file_rep, col_types = 'icd' )
    data <- bind_rows( data, data_rep )
}

# reorganize data for boxplots
# will appear in desired order (from kinship_methods)
data_list <- lapply( kinship_methods$code, function( x ) {
    # subset tibble to data from this kinship method only
    auc_x <- data[ data$kinship == x, ]
    # validate number of replicates
    stopifnot( nrow( auc_x ) == n_rep )
    stopifnot( all( auc_x$rep %in% 1 : n_rep ) )
    # just keep column of interest, only values of interest
    return( auc_x$herit )
})

# now make plot of data
dims <- fig_scale( 1.5 ) # w/h
fig_start(
    'herit',
    width = dims[1],
    height = dims[2],
    mar_b = 7
)
boxplot(
    data_list,
    names = NA, # individual labels will be plotted with rest
    xaxt = 'n',
    ylab = 'Heritability'
)
mtext( 'Kinship Estimate', side = 1, line = 6 )

abline( h = herit, lty = 2, col = 'blue' )
abline( h = herit_biased, lty = 2, col = 'red' )

# this hack makes it all more symmetrical!
kinship_methods$short[1] <- 'Popkin'
# reuse popkin-style labeling!  Great for hierarchical/factorial setup
# though popkin has defaults for these, the raw function doesn't have defaults!
popkin:::print_labels_multi(
             labs = cbind( kinship_methods$type, kinship_methods$short ),
             labs_cex = c(0.7, 0.7),
             labs_las = c(2, 0),
             labs_line = c(0.5, 4),
             labs_lwd = 1, # default
             labs_sep = c(FALSE, TRUE),
             labs_even = c(FALSE, TRUE),
             labs_ticks = FALSE, # default
             labs_text = TRUE, # default
             labs_col = 'black', # default
             # these align with barplot-specific changes
             xc_ind = 1L : nrow( kinship_methods ),
             doMat = FALSE
         )

fig_end()


