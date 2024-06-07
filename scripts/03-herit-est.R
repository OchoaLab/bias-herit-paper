# this script generates a few replicates of the same simulation with slight differences, to showcase the biases Zhuoran has identified.
# "mira2" forked with big removals/changes
# - removed all MVN tests, they no longer make sense (and combined tests with estimates from genotypes but MVN trait are even weirder)
# - separate genotype simulation/storage, precomputing and saving more than before!

library(optparse)    # for terminal options
library(tibble)      # to store data
library(readr)       # to write data
library(genbin)      # for running external binaries

# constants
# output data names (genotypes and trait)
name_phen <- 'data'

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# now move to data location
setwd( '../data/' )

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--name", type = "character", default = NA, 
                help = "Base name for genotype and phenotype simulation", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "replicate number", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (affects GCTA only)", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$name
rep <- opt$rep
threads <- opt$threads
m_causal <- opt$m_causal
herit <- opt$herit

# now move to data location
setwd( '../data/' )
# move into project
setwd( dir_out )

# dir must already exist since we need genotypes present!
setwd( paste0( 'rep-', rep ) )

# output path for gen/phen data
dir_out <- paste0(
  'mc', m_causal,
  '-h', herit
)

setwd( dir_out )


############
### GCTA ###
############

# wrapper for sapply "loop"
my_herit_lmm_gcta <- function ( name ) {
    name <- paste0( '../kinship/', name )
    herit <- gcta_reml(
        name, # kinship and output
        name_phen = name_phen,
        threads = threads
    )$herit
    # cleanup
    delete_files_gcta_hsq( name )
    delete_files_log( name )
    return( herit )
}

# compute and store directly into tibble
data <- tibble(
    rep = rep,
    kinship = kinship_methods$code,
    herit = herit,
    m_causal = m_causal,
    herit_est = sapply( kinship_methods$code, my_herit_lmm_gcta )
)

# save data
write_tsv( data, file = 'herit.txt.gz' )
