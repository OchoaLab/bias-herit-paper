# compare details of GWAS with true or biased kinship limits, and the two kinds of (noisy) estimates

library(optparse)    # for terminal options
library(readr)       # to write kinship matrix
library(tibble)      # to store data
library(genio)       # to write BED files for external software
library(popkin)      # to estimate kinship in LIGERA
library(popkinsuppl) # for PCA's kinship estimator

# switch to main scripts directory
setwd( '../../scripts/' )

# standard code for a complex trait and an admixed population
source('sim_geno_trait_k3.R')

# load new functions from external scripts
source('kinship_to_evd.R')
source('gas_pca_optim.R')
source('gas_lmm_gemma.R')

# get back to this subdirectory
setwd( '../bias/scripts/' )


# change this path to GEMMA to whatever it is in your computer
gemma_bin <- '/home/viiia/bin/gemma-0.98.1-linux-static'


############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
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
    make_option("--debug", action = "store_true", default = FALSE, 
                help = "debug mode (GEMMA is fully verbose)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
debug <- opt$debug

# output path for BED files
name_out <- paste0(
    'gas',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-g', generations
)

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
    generations = generations,
    herit = herit,
    verbose = TRUE,
    fst = fst
)
X <- obj$X
kinship <- obj$kinship
trait <- obj$trait
causal_indexes <- obj$causal_indexes
causal_coeffs <- obj$causal_coeffs

# write BED version for external code
plink_data <- write_plink(name_out, X, pheno = trait, verbose = FALSE)

###############
### KINSHIP ###
###############

# There are 4 kinship matrices to consider
# - `kinship`: true kinship matrix of simulation
# - `kinship_std_lim`: limit of biased "standard" estimator
# - `kinship_std`: biased (and noisy) "standard" estimate from genotypes (same as GEMMAs)
# - `kinship_popkin`: unbaised (but noisy) estimate from genotypes

# 1) true kinship (as given by simulation)

# 2) limit of biased "standard" estimator
# calculated with function kinship_std_limit from package popkinsuppl
kinship_std_lim <- kinship_std_limit( kinship )

# 3) compute biased "standard" kinship estimate
# calculated with function kinship_std from package popkinsuppl
kinship_std <- kinship_std( X )

# 4) popkin estimate
# need labels first
labs <- ceiling( ( 1 : n_ind ) / n_ind * 10 )
# actual popkin estimate
kinship_popkin <- popkin(X, labs)

#############
### GEMMA ###
#############

# code to write kinship tables for GEMMA
write_kinship_table <- function( kinship, file ) {
    kinship_tibble <- as_tibble( 2 * kinship, .name_repair = 'minimal')
    write_tsv(kinship_tibble, file, col_names = FALSE)
}

# test method with ideal matrix
message( "GEMMA (true kinship)" )
# write true kinship matrix in correct format
file_kinship_true <- paste0(name_out, '.kinship.txt')
write_kinship_table( kinship, file_kinship_true )
# run GEMMA
obj_gemma_true <- gas_lmm_gemma(gemma_bin, name_out, file_kinship_true, m_loci = m_loci, debug = debug)
## obj_gemma_true$pvals
## obj_gemma_true$beta_hat

# test method with popkin estimate
message( "GEMMA (popkin kinship)" )
# write true kinship matrix in correct format
file_kinship_popkin <- paste0(name_out, '.kinship_popkin.txt')
write_kinship_table( kinship_popkin, file_kinship_popkin )
# run GEMMA
obj_gemma_popkin <- gas_lmm_gemma(gemma_bin, name_out, file_kinship_popkin, m_loci = m_loci, debug = debug)
## obj_gemma_popkin$pvals
## obj_gemma_popkin$beta_hat

# GEMMA with limit of biased "standard" kinship estimate
message( "GEMMA (std lim kinship)" )
# write true kinship matrix in correct format
file_kinship_std_lim <- paste0(name_out, '.kinship_std_lim.txt')
write_kinship_table( kinship_std_lim, file_kinship_std_lim )
# run GEMMA
obj_gemma_std_lim <- gas_lmm_gemma(gemma_bin, name_out, file_kinship_std_lim, m_loci = m_loci, debug = debug)
## obj_gemma_std_lim$pvals
## obj_gemma_std_lim$beta_hat

# GEMMA with biased "standard" kinship estimate
message( "GEMMA (std kinship)" )
## # obtain standard kinship estimate using GEMMA itself
## file_gemma_kin <- gas_lmm_gemma_kin(gemma_bin, name_out, debug = debug)$file
# write true kinship matrix in correct format
file_kinship_std <- paste0(name_out, '.kinship_std.txt')
write_kinship_table( kinship_std, file_kinship_std )
# run GEMMA
obj_gemma_std <- gas_lmm_gemma(gemma_bin, name_out, file_kinship_std, m_loci = m_loci, debug = debug)
## obj_gemma_std$pvals
## obj_gemma_std$beta_hat

# cleanup: do after all GEMMA runs are finished
# remove gemma's kinship estimate files
invisible( file.remove( file_kinship_true ) )
invisible( file.remove( file_kinship_popkin ) )
invisible( file.remove( file_kinship_std ) )
invisible( file.remove( file_kinship_std_lim ) )
# delete the temporary BED files now that we are done
delete_files_plink(name_out)

###########
### PCA ###
###########

# indexes for PCs
# true version and popkin uses all
indexes_true <- 1 : k_subpops
# std versions uses K-1 (theory suggests this makes more sense)
indexes_std <- 1 : ( k_subpops - 1 )

message( "PCA (true kinship)" )
# oracle version uses true kinship matrix
eigenvectors <- kinship_to_evd( kinship )
obj_pca_true <- gas_pca_optim(X, trait, eigenvectors[, indexes_true ])
## obj_pca_true$pvals
## obj_pca_true$beta_hat

message( "PCA (popkin kinship)" )
eigenvectors_popkin <- kinship_to_evd( kinship_popkin )
obj_pca_popkin <- gas_pca_optim(X, trait, eigenvectors_popkin[, indexes_true ])
## obj_pca_popkin$pvals
## obj_pca_popkin$beta_hat

message( "PCA (std lim kinship)" )
eigenvectors_std_lim <- kinship_to_evd( kinship_std_lim )
obj_pca_std_lim <- gas_pca_optim(X, trait, eigenvectors_std_lim[, indexes_std ])
## obj_pca_std_lim$pvals
## obj_pca_std_lim$beta_hat

message( "PCA (std kinship)" )
eigenvectors_std <- kinship_to_evd( kinship_std )
obj_pca_std <- gas_pca_optim(X, trait, eigenvectors_std[, indexes_std ])
## obj_pca_std$pvals
## obj_pca_std$beta_hat


