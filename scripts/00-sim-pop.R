# simulates overall population structure (info shared across reps, which are created later)

library(optparse)    # for terminal options
library(simgenphen)  # wrapper for various simulation packages
library(genio)       # to write files for external software
library(popkinsuppl) # for kinship_std estimator

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 10, 
                help = "number of replicates", metavar = "int"),
    make_option(c("-n", "--n_ind"), type = "integer", default = 2000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 3, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.3, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_rep <- opt$rep
n_ind <- opt$n_ind
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
G <- opt$generations

# output path for BED files
dir_out <- paste0(
    'sim',
    '-n', n_ind,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-g', G
)

# now move to data location
setwd( '../data/' )

# dir shouldn't exist because this is run only once (do not overwrite!)
if ( dir.exists( dir_out ) )
    stop( 'Data exists, will not ovewrite: ', dir_out )

# create and move in
dir.create( dir_out )
setwd( dir_out )

# the population structure is the same across replicates
obj <- sim_pop(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    fst = fst,
    G = G
)
admix_proportions_1 <- obj$admix_proportions_1
inbr_subpops <- obj$inbr_subpops
famG <- obj$fam
ids <- obj$ids
kinship <- obj$kinship

# save simulation details to file
# (kinship is stored as GRM and not needed for sims anyway)
save( admix_proportions_1, inbr_subpops, famG, ids, file = 'simpop.RData' )

# put all kinship matrices in a subdirectory
dir.create( 'kinship' )
setwd( 'kinship' )

# true kinship (as given by simulation)
write_grm( 'true', 2 * kinship )

# limit of biased "standard" estimator (from package popkinsuppl)
write_grm( 'std_rom_lim', 2 * kinship_std_limit( kinship ) )


