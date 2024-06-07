# simulate phenotypes in reps

library(optparse)    # for terminal options
library(simgenphen)  # wrapper for various simulation packages
library(genio)       # to write files for external software

# constants
# output data names (genotypes and phenotypes)
name_out <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--name", type = "character", default = NA, 
                help = "Base name for population structure simulation", metavar = "character"),
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

# now move to data location
setwd( '../data/' )
# move into project
setwd( opt$name )




rep_dir <- paste0( 'rep-', rep )

setwd( rep_dir )

load( 'p_anc.Rdata' )

obj<-read_plink( name_out )

X<-obj$X
fam<-obj$fam

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

# and the proper genetic traits
obj <- sim_trait_env(
    X = X,
    m_causal = m_causal,
    herit = herit,
    p_anc = p_anc
)

fam$pheno<-obj$trait

# write to phen for later use
write_phen( name_out, fam )
