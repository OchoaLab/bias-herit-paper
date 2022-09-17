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
# for non-truth, should technically use every rep, though here we'll assume this is a stable value
mean_kinship_popkin_rom <- mean( read_grm( 'rep-1/kinship/popkin_rom' )$kinship / 2 )
mean_kinship_popkin_mor <- mean( read_grm( 'rep-1/kinship/popkin_mor' )$kinship / 2 )
# actually predict biased herit value from formula
herit_biased_true <- herit * ( 1 - mean_kinship_true ) / ( 1 - ( mean_kinship_true * herit ) )
herit_biased_popkin_rom <- herit * ( 1 - mean_kinship_popkin_rom ) / ( 1 - ( mean_kinship_popkin_rom * herit ) )
herit_biased_popkin_mor <- herit * ( 1 - mean_kinship_popkin_mor ) / ( 1 - ( mean_kinship_popkin_mor * herit ) )

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
    ylab = 'Heritability estimate'
)
mtext( 'Kinship estimate', side = 1, line = 6 )

# lines and legend
abline( h = herit, lty = 2, col = 'blue' )
abline( h = herit_biased_true, lty = 2, col = 'red' )
abline( h = herit_biased_popkin_rom, lty = 2, col = 'pink' )
abline( h = herit_biased_popkin_mor, lty = 2, col = 'green' )
legend(
    'bottomleft',
    c('Truth', 'Bias ROM lim.', 'Bias ROM est.', 'Bias MOR est.'),
    lty = 2,
    col = c('blue', 'red', 'pink', 'green'),
    cex = 0.7
)

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


# make smaller "grant proposal" version with fewer data
# define what to keep
codes_small <- c('true', 'popkin_rom', 'popkin_mor', 'std_mor')
indexes <- kinship_methods$code %in% codes_small
data_list_small <- data_list[ indexes ]
kinship_methods_small <- kinship_methods[ indexes, ]
# edit names some more (space is precious)
kinship_methods_small$nice[ kinship_methods_small$code == 'true' ] <- 'True'
kinship_methods_small$nice <- sub( ' est.', '', kinship_methods_small$nice )
# here we'll go with normal names under boxplots
names( data_list_small ) <- kinship_methods_small$nice

fig_start(
    'herit-small',
    width = 2.5, # way smaller than full fig, for grant
    height = 4,
    mar_b = 5.5
)
# to control labels underneath plot
par_orig <- par( las = 3, cex.axis = 0.7 )
# actual plot
boxplot(
    data_list_small,
    ylab = 'Heritability estimate'
)
# restore `las` and anything else that might have gotten messed up
par( par_orig )
mtext( 'Kinship estimate', side = 1, line = 4.5 )

# lines and legend
abline( h = herit, lty = 2, col = 'blue' )
abline( h = herit_biased_popkin_mor, lty = 2, col = 'red' )
legend(
    'bottomleft',
    c('Truth', 'Predicted bias'),
    lty = 2,
    col = c('blue', 'red'),
    cex = 0.7
)

fig_end()
