# compare details of GWAS with true or biased kinship limits, and their noisy estimates
library(optparse)    # for terminal options
library(readr)       # to write kinship matrix
library(genio)       # to write BED files for external software
library(popkin)      # to estimate kinship in RGLS
library(popkinsuppl) # for PCA's kinship estimator

# switch to main scripts directory
setwd( '../../scripts/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/gas-rgls/scripts' )

# standard code for a complex trait and an admixed population
source('sim_geno_trait_k3.R')

# load new functions from external scripts
source('kinship_to_evd.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')
source('kinship_gcta_limit.R')

# place outputs in "data" (in a subdirectory depending on params)
setwd( '../bias/data/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/true-vs-biased-kinship-gwas/results' )

# change this path to GCTA to whatever it is in your computer
gcta_bin <- '/home/viiia/bin/gcta_1.93.2beta/gcta64'
#gcta_bin <- 'D:/3.Duke/research/alex_ochoa/software/gcta_1.93.2beta_win/bin/gcta64'

# a name for temporary BED/etc data, under project dir
name_out <- 'data'

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 10000, 
                help = "number of loci", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 3, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.3, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal", type = "integer", default = 100, 
                help = "num causal loci", metavar = "int"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (affects GCTA only)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
G <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
fes <- opt$fes
threads <- opt$threads


# output path for BED files and all results files
dir_out <- paste0(
    'sim-admix',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-g', G,
    if ( fes ) '-fes' else '-rc'
)


# place data in a new figure specific to this simulation
if ( !dir.exists( dir_out ) )
    dir.create( dir_out )
setwd( dir_out )


############
### SIMS ###
############

# simulate genotypes and trait as usual
obj <- sim_geno_trait_k3(
    n_ind = n_ind,
    m_loci = m_loci,
    m_causal = m_causal,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    G = G,
    herit = herit,
    fes = fes,
    verbose = TRUE,
    fst = fst
)
X <- obj$X
kinship_true <- obj$kinship
trait <- obj$trait
causal_indexes <- obj$causal_indexes
causal_coeffs <- obj$causal_coeffs

# save causal trait data for AUC evaluations later
save(
    trait,
    causal_indexes,
    causal_coeffs,
    file = 'simtrait.RData'
)

# write BED version for external code
plink_data <- write_plink(name_out, X, pheno = trait, verbose = FALSE)
# write phenotype file
write_phen(name_out, plink_data$fam, verbose = FALSE)


###############
### KINSHIP ###
###############

# There are all of these kinship matrices to consider
# - `kinship_true`: true kinship matrix of simulation
# - `kinship_std_rom`: biased (and noisy) ROM version of "standard" estimate from genotypes
# - `kinship_std_rom_lim`: limit of biased ROM version of "standard" estimator
# - `kinship_std_mor`: biased (and noisy) MOR version of "standard" estimate from genotypes (same as GEMMAs)
# - `kinship_popkin`: unbaised (but noisy) estimate from genotypes
# - `kinship_wg`: biased (and noisy) WG estimate from genotypes
# - `kinship_wg_lim`: limit of biased WG estimator
# - `kinship_gcta`: biased (and noisy) GCTA estimate from genotypes
# - `kinship_gcta_lim`: limit of biased GCTA estimator

# 1) true kinship (as given by simulation)

# 2) limit of biased "standard" estimator
# calculated with function kinship_std_limit from package popkinsuppl
kinship_std_rom_lim <- kinship_std_limit( kinship_true )

# 3) compute biased "standard" kinship estimate
# calculated with function kinship_std from package popkinsuppl
kinship_std_rom <- kinship_std( X )
kinship_std_mor <- kinship_std( X, mean_of_ratios = TRUE )

# 4) popkin estimate
# need labels first
labs <- ceiling( ( 1 : n_ind ) / n_ind * 10 )
# actual popkin estimate
kinship_popkin <- popkin(X, labs)

# 5) WG
# limit of estimator
kinship_wg_lim <- kinship_wg_limit( kinship_true )
# estimate
kinship_wg <- kinship_wg_limit( kinship_popkin )

# 6) GCTA
# actual estimate is given directly by GCTA, below
# limit of estimator
kinship_gcta_lim <- kinship_gcta_limit( kinship_true )
# estimate kinship with GCTA's method
gas_lmm_gcta_kin(gcta_bin, name_out, name_out = 'kinship_gcta')
unlink( 'kinship_gcta.log' ) # unneeded log file
# read GCTA kinship into R (need for PCA version)
obj <- read_grm( 'kinship_gcta', verbose = FALSE )
kinship_gcta <- obj$kinship / 2
M <- obj$M # same for all methods (control that aspect)

############
### GCTA ###
############

# NOTE: for making plots/etc, will preserve all GRMs!

message( "lmm_gcta" )
# NOTE: GCTA's kinship estimate was produced earlier already!
# GWAS using its own kinship matrix
obj <- gas_lmm_gcta(
    gcta_bin,
    name_out,
    name_grm = 'kinship_gcta',
    m_loci = m_loci,
    threads = threads
)
# cleanup
delete_files_gcta(name_out) # GAS table
# initialize these data frames with output data
pvals <- data.frame( lmm_gcta = obj$pvals )
betas <- data.frame( lmm_gcta = obj$beta_hat )

message( "pca_gcta" )
# indexes for PCs
indexes <- 1 : k_subpops
eigenvectors <- kinship_to_evd( kinship_gcta )
obj <- gas_pca_optim(X, trait, eigenvectors[, indexes ])
# add to data frames
pvals$pca_gcta <- obj$pvals
betas$pca_gcta <- obj$beta_hat

# automated version for all other cases (only GCTA with default kinship matrix makes more sense differently)
do_all <- function(kinship) {
    # use actual variable name as file output (all is consistent for us)
    name_grm <- deparse( substitute( kinship ) )
    name_method <- sub( '^kinship_', '', name_grm ) # don't need "kinship_" prefix in outputs
    name_lmm_method <- paste0( 'lmm_', name_method ) # for LMM outputs
    name_pca_method <- paste0( 'pca_', name_method ) # for PCA outputs

    # LMM
    message( name_lmm_method )
    # use shared globals and typical processing
    write_grm(
        name_grm,
        kinship = 2 * kinship,
        M = M,
        fam = plink_data$fam,
        verbose = FALSE
    )
    # GWAS
    obj <- gas_lmm_gcta(
        gcta_bin,
        name_out,
        name_grm = name_grm,
        m_loci = m_loci,
        threads = threads
    )
    # cleanup
    delete_files_gcta( name_out ) # GAS table
    # add to data frames
    pvals[[ name_lmm_method ]] <- obj$pvals
    betas[[ name_lmm_method ]] <- obj$beta_hat
    
    # PCA
    message( name_pca_method )
    eigenvectors <- kinship_to_evd( kinship )
    obj <- gas_pca_optim(X, trait, eigenvectors[, indexes ])
    # add to data frames
    pvals[[ name_pca_method ]] <- obj$pvals
    betas[[ name_pca_method ]] <- obj$beta_hat

    # make sure globals are edited
    pvals <<- pvals
    betas <<- betas
}

do_all( kinship_gcta_lim )
do_all( kinship_true )
do_all( kinship_popkin )
do_all( kinship_std_rom )
do_all( kinship_std_rom_lim )
do_all( kinship_std_mor )
do_all( kinship_wg )
do_all( kinship_wg_lim )

# final cleanup: do after all GCTA runs are finished
# delete the temporary BED files now that we are done
delete_files_plink(name_out)
delete_files_phen(name_out)

# save data frames!
write_tsv( pvals, 'pvals.txt' )
write_tsv( betas, 'betas.txt' )
