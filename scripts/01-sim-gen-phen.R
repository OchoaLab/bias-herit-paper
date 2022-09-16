# simulate genotypes and phenotypes in reps

library(optparse)    # for terminal options
library(simgenphen)  # wrapper for various simulation packages
library(genio)       # to write files for external software

# constants
# output data names (genotypes and trait)
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
    make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
                help = "number of loci", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
rep <- opt$rep
m_loci <- opt$m_loci
m_causal <- opt$m_causal
herit <- opt$herit

# now move to data location
setwd( '../data/' )
# move into project
setwd( opt$name )
# load population structure variables
load( 'simpop.RData' )
# loads: admix_proportions_1, inbr_subpops, famG, ids

# output path for gen/phen data
dir_out <- paste0(
    'm', m_loci,
    '-mc', m_causal,
    '-h', herit
)

# may exist (if adding reps), but create otherwise
if ( !dir.exists( dir_out ) )
    dir.create( dir_out )
setwd( dir_out )

rep_dir <- paste0( 'rep-', rep )
# don't overwrite existing data
if ( dir.exists( rep_dir ) )
    stop( 'Rep exists: ', rep_dir )

# else create and move into dir
dir.create( rep_dir )
setwd( rep_dir )

################
### SIM GENO ###
################

# simulate new genotypes for each replicate
obj <- sim_geno(
    admix_proportions_1 = admix_proportions_1,
    inbr_subpops = inbr_subpops,
    fam = famG,
    ids = ids,
    m_loci = m_loci
)
X <- obj$X
p_anc <- obj$p_anc

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
# causal_indexes <- obj$causal_indexes # unused
# causal_coeffs <- obj$causal_coeffs # unused

# write to bed/bim/fam/phen for later use
write_plink( name_out, X, pheno = obj$trait, write_phen = TRUE )
