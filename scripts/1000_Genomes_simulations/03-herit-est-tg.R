
library(optparse)    # for terminal options
library(tibble)      # to store data
library(readr)       # to write data
library(genbin)      # for running external binaries
library(tidyverse)

#devtools::install_github('ochoalab/genbin')


############
### ARGV ###
############

# define options
option_list = list(
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

rep <- opt$rep
threads <- opt$threads
m_causal <- opt$m_causal
herit <- opt$herit

# now move to data location
setwd( 'traits_path')


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

data <- NULL
data_row<- NULL

for ( grm in c('tgp_qc_std_mor','tgp_qc_std_rom','tgp_qc_popkin_rom')) {

for (grm_maf in c('maf000','maf001','maf005')) {

for (m_pheno in 1:3) {
# wrapper for sapply "loop"
 herit_est <- gcta_reml(
        'traits', 
        name_grm = paste0( '../../../', grm_maf, '/', grm ),
        threads = threads,
        m_pheno = m_pheno
    )$herit
    # cleanup
delete_files_gcta_hsq( 'traits' )
delete_files_log( 'traits' )

data_row <- tibble(
  rep = rep,
  kinship = grm,
  herit = herit,
  m_causal = m_causal,
  grm_maf = grm_maf,
  m_pheno = m_pheno,
  herit_est = herit_est
)

data <- bind_rows( data, data_row )

}
}
}

# save data
write_tsv( data, file = 'herit.txt.gz' )
