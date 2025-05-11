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
                help = "Total number of replicates", metavar = "int"),
    make_option("--n_h", type = "integer", default = NA, 
                help = "Total number of heritability", metavar = "int"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$name
n_rep <- opt$n_rep
n_h<-opt$n_h
m_causal <- opt$m_causal


if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

if ( is.na( n_h ) )
  stop( 'Option `--n_h` is required!' )

# go where the data is
setwd( '../data/' )
setwd( dir_out )


## create folder list
h_list<-c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

dir_out_list<-NULL
# output path for gen/phen data

i=1
for (herit in h_list ) {
  dir_out_list[i] <- paste0(
    'mc', m_causal,
    '-h', herit
  )
  i=i+1
  
}

# tibble to grow
data <- NULL

# load pre-calculated herit estimates
for ( rep in 1 : n_rep ) {
  
  for (h in 1: n_h) {
  
    file_rep <- paste0( 'rep-', rep,'/',dir_out_list[h], '/herit.txt.gz' )
    data_rep <- read_tsv( file_rep, col_types = 'icd' )
    data <- bind_rows( data, data_rep )
    
  }
}

# save into combined table
write_tsv( data, 'herit.txt.gz' )




