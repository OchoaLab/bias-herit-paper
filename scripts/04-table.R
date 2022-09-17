# gathers small tables from each rep into a single big one (for easy commits and downstream analysis)

library(optparse) # for terminal options
library(readr)    # to read tables
library(dplyr)    # for bind_rows

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

# go where the data is
setwd( '../data/' )
setwd( dir_out )

# tibble to grow
data <- NULL

# load pre-calculated herit estimates
for ( rep in 1 : n_rep ) {
    file_rep <- paste0( 'rep-', rep, '/herit.txt.gz' )
    data_rep <- read_tsv( file_rep, col_types = 'icd' )
    data <- bind_rows( data, data_rep )
}

# save into combined table
write_tsv( data, 'herit.txt.gz' )
