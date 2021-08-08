# this script generates a few replicates of the same simulation with slight differences, to showcase the biases Zhuoran has identified.

library(optparse)    # for terminal options
library(tibble)      # to store data
library(readr)       # to write data
library(genio)       # to write files for external software
library(popkin)      # for popkin
library(popkinsuppl) # for kinship_std estimator

# switch to main scripts directory
dir_orig <- getwd() # remember where we are now
setwd( '../../scripts/' )
source('sim_geno_trait_k3.R')
source('gas_lmm_gcta.R')
source('kinship_gcta_limit.R')
#source('herit_lmmlite.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# now move to data location
setwd( '../data/' )

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 10, 
                help = "number of replicates", metavar = "int"),
    make_option(c("-n", "--n_ind"), type = "integer", default = 2000, 
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
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (affects GCTA only)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
rep <- opt$rep
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
G <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
threads <- opt$threads

# output path for BED files
name_out <- paste0(
    'herit-estimates',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-g', G,
    '-r', rep
)

###############
### KINSHIP ###
###############

# the population structure is the same across replicates
obj <- sim_pop(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    fst = fst,
    G = G,
    verbose = TRUE
)
admix_proportions <- obj$admix_proportions
inbr_subpops <- obj$inbr_subpops
famG <- obj$fam
ids <- obj$ids
kinship <- obj$kinship

# There are 6 kinship matrices to consider
# 1- `kinship`: true kinship matrix of simulation
# 2- `kinship_std_lim`: limit of biased "standard" estimator
# 3- `kinship_gcta_lim`: limit of the GCTA estimator
# 4- `kinship_popkin`: unbiased (but noisy) estimate from genotypes
# 5- `kinship_std`: biased (and noisy) "standard" estimate from genotypes (same as GEMMAs)
# 6- `kinship_gcta`: GCTA estimator (equation is most similar to `std`)

# 1) true kinship (as given by simulation)

# 2) limit of biased "standard" estimator
# calculated with function kinship_std_limit from package popkinsuppl
kinship_std_lim <- kinship_std_limit( kinship )

# 3) gcta limit estimate
kinship_gcta_lim <- kinship_gcta_limit( kinship )

# SAVE as GRM files

# matrix of sample sizes is complete matrix in all these cases (no missingness)
M <- matrix( m_loci, nrow = n_ind, ncol = n_ind )
# also need a dummy fam file for writing the GRM
fam <- make_fam( n = n_ind )

# to avoid mistakes, define a smaller wrapper
my_write_grm <- function( name_out, kinship ) {
    write_grm(
        name_out,
        kinship = 2 * kinship,
        M = M,
        fam = fam,
        verbose = TRUE
    )
}



# write the kinship matrix to a file
name_kinship_true <- paste0(name_out, '.kinship_true')
my_write_grm( name_kinship_true, kinship )

name_kinship_std_lim <- paste0(name_out, '.kinship_std_lim')
my_write_grm( name_kinship_std_lim, kinship_std_lim )

name_kinship_gcta_lim <- paste0(name_out, '.kinship_gcta_lim')
my_write_grm( name_kinship_gcta_lim, kinship_gcta_lim )

#################
### SIM TRAIT ###
#################

message( 'sim_trait_mvn...' )

# draw all MVN traits at once (most efficient)
traits_mvn <- sim_trait_mvn(
    rep = rep,
    kinship = kinship,
    herit = herit
)

###################
### GENOME REPS ###
###################

# store the heritability estimates in these vectors
# first using trait from genotypes
herit_truX <- vector( 'numeric', rep )
herit_popX <- vector( 'numeric', rep )
herit_stdX <- vector( 'numeric', rep )
herit_slmX <- vector( 'numeric', rep )
herit_glmX <- vector( 'numeric', rep )
herit_gctX <- vector( 'numeric', rep )
# then using trait from MVN model
herit_truN <- vector( 'numeric', rep )
herit_popN <- vector( 'numeric', rep )
herit_stdN <- vector( 'numeric', rep )
herit_slmN <- vector( 'numeric', rep )
herit_glmN <- vector( 'numeric', rep )
herit_gctN <- vector( 'numeric', rep )

for ( rep_i in 1 : rep ) {
    
    ################
    ### SIM GENO ###
    ################

    # simulate new genotypes for each replicate
    obj <- sim_geno(
        admix_proportions = admix_proportions,
        inbr_subpops = inbr_subpops,
        G = G,
        fam = famG,
        ids = ids,
        m_loci = m_loci,
        verbose = TRUE
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
        p_anc = p_anc,
        verbose = TRUE
    )
    # since we have two traits, differentiate in names
    trait_X <- obj$trait
    causal_indexes <- obj$causal_indexes
    causal_coeffs <- obj$causal_coeffs
    
    # write phenotype file (this one is the true genetic trait)
    fam$pheno <- trait_X
    name_pheno_X <- paste0(name_out, '.pheno_X')
    write_phen(name_pheno_X, fam)

    # and write MVN trait too
    fam$pheno <- traits_mvn[ rep_i, ]
    name_pheno_N <- paste0(name_out, '.pheno_N')
    write_phen(name_pheno_N, fam)
    
    # write BED version for external code
    write_plink(name_out, X)
    
    ###############
    ### KINSHIP ###
    ###############

    # estimates that get redone at every replicate

    # 4) popkin estimate
    # need labels first
    labs <- ceiling( ( 1 : n_ind ) / n_ind * 10 )
    # actual popkin estimate
    message( 'popkin' )
    kinship_popkin <- popkin(X, labs)

    # 5) compute biased "standard" kinship estimate
    # calculated with function kinship_std from package popkinsuppl
    message( 'kinship_std' )
    kinship_std <- kinship_std( X )

    # SAVE as GRM files
    name_kinship_popkin <- paste0(name_out, '.kinship_popkin')
    my_write_grm( name_kinship_popkin, kinship_popkin )
    
    name_kinship_std <- paste0(name_out, '.kinship_std')
    my_write_grm( name_kinship_std, kinship_std )

    # 6) GCTA is created by its own software
    name_kinship_gcta <- paste0(name_out, '.kinship_gcta')
    gas_lmm_gcta_kin(gcta_bin, name = name_out, name_out = name_kinship_gcta, threads = threads)
    # cleanup
    delete_files_log(name_kinship_gcta)
    
    ## #########################
    ## ### LMMLITE PREROTATE ###
    ## #########################

    ## message( 'herit_lmmlite_prerotate...' )

    ## # saves time factoring an eigendecomposition for LMMLITE
    ## prerot_kinship <- herit_lmmlite_prerotate( kinship )
    ## prerot_kinship_std <- herit_lmmlite_prerotate( kinship_std )
    ## prerot_kinship_std_lim <- herit_lmmlite_prerotate( kinship_std_lim )
    ## prerot_kinship_popkin <- herit_lmmlite_prerotate( kinship_popkin )

    ## ###############
    ## ### lmmlite ###
    ## ###############

    ## for ( i in 1 : rep ) {
    ##     trait <- traits[ i, ]
    
    ##     message( "lmmlite (true kinship)" )
    ##     herit_all[i,1] <- herit_lmmlite(trait, prerot_kinship, prerotated = TRUE)

    ##     message( "lmmlite (popkin kinship)" )
    ##     herit_all[i,2] <- herit_lmmlite(trait, prerot_kinship_popkin, prerotated = TRUE)

    ##     #lmmlite with limit of biased "standard" kinship estimate
    ##     message( "lmmlite (std lim kinship)" )
    ##     herit_all[i,3] <- herit_lmmlite(trait, prerot_kinship_std_lim, prerotated = TRUE)

    ##     # lmmlite with biased "standard" kinship estimate
    ##     message( "lmmlite (std kinship)" )
    ##     herit_all[i,4] <- herit_lmmlite(trait, prerot_kinship_std, prerotated = TRUE)
    ## }

    ############
    ### GCTA ###
    ############

    # get estimates
    # another wrapper to avoid mistakes in these loops
    my_herit_lmm_gcta <- function ( name_pheno, name_kinship ) {
        herit_lmm_gcta(gcta_bin, name_pheno, name_grm = name_kinship, threads = threads)$herit
    }

    # first using trait from genotypes
    
    message('herit_truX...')
    herit_truX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_true)
    
    message('herit_popX...')
    herit_popX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_popkin)
    
    message('herit_stdX...')
    herit_stdX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_std)
    
    message('herit_slmX...')
    herit_slmX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_std_lim)

    message('herit_glmX...')
    herit_glmX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_gcta_lim)
    
    message('herit_gctX...')
    herit_gctX[ rep_i ] <- my_herit_lmm_gcta(name_pheno_X, name_kinship_gcta)
    
    # then using trait from MVN model
    
    message('herit_truN...')
    herit_truN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_true)
    
    message('herit_popN...')
    herit_popN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_popkin)
    
    message('herit_stdN...')
    herit_stdN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_std)
    
    message('herit_slmN...')
    herit_slmN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_std_lim)

    message('herit_glmN...')
    herit_glmN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_gcta_lim)
    
    message('herit_gctN...')
    herit_gctN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_gcta)
    
    ###############
    ### CLEANUP ###
    ###############
    
    # these get deleted after every iteration
    delete_files_plink(name_out)
    delete_files_phen(name_pheno_X)
    delete_files_phen(name_pheno_N)
    delete_files_grm(name_kinship_popkin)
    delete_files_grm(name_kinship_std)
    delete_files_grm(name_kinship_gcta)
    # herit output and log files (only two, they get overwritten for every GRM tested)
    delete_files_gcta(name_pheno_X, herit = TRUE)
    delete_files_gcta(name_pheno_N, herit = TRUE)
}

######################
### SAVE ESTIMATES ###
######################

# gather vectors into a tibble
data <- tibble(
    truX = herit_truX,
    popX = herit_popX,
    slmX = herit_slmX,
    stdX = herit_stdX,
    glmX = herit_glmX,
    gctX = herit_gctX,
    truN = herit_truN,
    popN = herit_popN,
    slmN = herit_slmN,
    stdN = herit_stdN,
    glmN = herit_glmN,
    gctN = herit_gctN
)

# save data
write_tsv(
    data,
    file = paste0( name_out, '.txt' )
)

###############
### CLEANUP ###
###############

# delete these shared files only after all replicates are done
delete_files_grm(name_kinship_true)
delete_files_grm(name_kinship_std_lim)
delete_files_grm(name_kinship_gcta_lim)
