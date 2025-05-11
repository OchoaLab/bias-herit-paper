# simulate phenotypes in reps

library(optparse)    # for terminal options
library(simtrait)  # wrapper for various simulation packages
library(genio)       # to write files for external software
library(BEDMatrix)
library(tidyverse)

#devtools::install_github('ochoalab/simtrait')

# constants
# output data names (genotypes and phenotypes)
genotype_path <- 'genotype_path'


############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "replicate number", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
rep <- opt$rep
m_causal <- opt$m_causal
herit <- opt$herit

# move into project

setwd( 'traits_path')

rep_dir <- paste0( 'rep-', rep )

if ( !dir.exists( rep_dir ) )
  dir.create( rep_dir )

setwd( rep_dir )


X<-BEDMatrix( genotype_path )

fam<-read_fam(genotype_path)
fam<-fam[,1:2]

# output path for gen/phen data
dir_out <- paste0(
  'mc', m_causal,
  '-h', herit
)

# may exist (if adding reps), but create otherwise
if ( !dir.exists( dir_out ) )
  dir.create( dir_out )
setwd( dir_out )



#################
### SIM TRAIT ###
#################

for (maf_cut_value in c(0,0.01,0.05) ) {

# and the proper genetic traits
obj <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    kinship = switch(as.character(maf_cut_value),
                     "0" = 0.133,
                     "0.01" = 0.134,
                     "0.05"  = 0.1,
                     NA),
    maf_cut= maf_cut_value,
    fes = TRUE, fes_kinship_method = 'mle'
)

fam<-bind_cols(fam,obj$trait)

}

# write to phen for later use
write_tsv( fam,  'traits.phen', col_names=FALSE )


