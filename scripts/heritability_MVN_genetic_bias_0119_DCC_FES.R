# this script generates a few replicates of the same simulation with slight differences, to showcase the biases Zhuoran has identified.

library(optparse)    # for terminal options
library(tibble)      # to store data
library(readr)       # to write data
library(genio)       # to write files for external software
library(popkin)      # for popkin
library(popkinsuppl) # for kinship_std estimator


# switch to main scripts directory
dir_orig <- getwd() # remember where we are now
#setwd( '../../scripts/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/gas-rgls/scripts' )
#setwd('/home/dave/Research_Alex/1.reverse_regression/coding/mycode/herit_1224_DCC/DCC_version_0103')
setwd('./scripts')


source('sim_geno_trait_k3.R')
source('gas_lmm_gcta.R')
source('kinship_gcta_limit.R')
#source('herit_lmmlite.R')
#source('paths.R')
#setwd('/home/dave/Research_Alex/1.reverse_regression/coding/mycode/herit_1224_DCC/DCC_version_0103')
setwd( dir_orig ) # go back to where we were

# now move to data location
setwd( './results_mvn_genetic_bias3/')


# change this path to GCTA to whatever it is in your computer
#gcta_bin <- '/home/viiia/bin/gcta_1.93.2beta/gcta64'
#gcta_bin <- 'D:/3.Duke/research/alex_ochoa/software/gcta_1.93.2beta_win/bin/gcta64'
#gcta_bin <- '/home/dave/Research_Alex/software/gcta_1.93.2beta/gcta64'
#gcta_bin <- './software/gcta_1.93.2beta/gcta64'
gcta_bin <- '/hpc/home/zh105/Alex_project/gcta_1.93.2beta/gcta64'

#p_anc_input <- NULL
p_anc_input <- 0.5

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "number of replicates", metavar = "int"),
    make_option(c("-n", "--n_ind"), type = "integer", default = 2000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 200000, 
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
    make_option("--m_causal", type = "integer", default = 200000, 
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

# There are 1 kinship matrices to consider
# 1- `kinship`: true kinship matrix of simulation

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

loci_test <- c(1000,10000,100000)
#loci_test <- c(100000)

# store the heritability estimates in these vectors
# first using trait from genotypes
#herit_truX <- vector( 'numeric', rep )
herit_truX1 <- matrix( NA, ncol = length(loci_test), nrow = rep )
herit_truX2 <- matrix( NA, ncol = length(loci_test), nrow = rep )

# then using trait from MVN model
herit_truN <- vector( 'numeric', rep )



##############################################
### Estimate heritability based on MVN trait##
##############################################

for ( rep_i in 1 : rep ) {
    
    # and write MVN trait too
    fam$pheno <- traits_mvn[ rep_i, ]
    name_pheno_N <- paste0(name_out, '.pheno_N')
    write_phen(name_pheno_N, fam)
    
    ############
    ### GCTA ###
    ############
    
    # get estimates
    # another wrapper to avoid mistakes in these loops
    my_herit_lmm_gcta <- function ( name_pheno, name_kinship ) {
        herit_lmm_gcta(gcta_bin, name_pheno, name_grm = name_kinship, threads = threads)$herit
    }
    
    # then using trait from MVN model
    message('herit_truN...')
    herit_truN[ rep_i ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_true)
    
    
    ###############
    ### CLEANUP ###
    ###############
    
    # these get deleted after every iteration
    delete_files_phen(name_pheno_N)
    
    
    # herit output and log files (only two, they get overwritten for every GRM tested)
    delete_files_gcta(name_pheno_N, herit = TRUE)
}

##################################################
### Estimate heritability based on genetic trait##
##################################################

i_loci <-1

for ( m_loci in loci_test  ) {
    
    m_causal <- m_loci

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
        verbose = TRUE,
        p_anc_input = p_anc_input
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
        verbose = TRUE,
        fes = TRUE
    )
    # since we have two traits, differentiate in names
    trait_X1 <- obj$trait

    # and the proper genetic traits
    obj2 <- sim_trait_env(
        X = X,
        m_causal = m_causal,
        herit = herit,
        p_anc = p_anc,
        verbose = TRUE,
        fes = FALSE
    )
    # since we have two traits, differentiate in names
    trait_X2 <- obj2$trait
    
    
    # write phenotype file (this one is the true genetic FES)
    fam$pheno <- trait_X1
    name_pheno_X1 <- paste0(name_out, '.pheno_X1')
    write_phen(name_pheno_X1, fam)
    
    # write phenotype file (this one is the true genetic RC)
    fam$pheno <- trait_X2
    name_pheno_X2 <- paste0(name_out, '.pheno_X2')
    write_phen(name_pheno_X2, fam)

    # write BED version for external code
    #write_plink(name_out, X)
    

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
    herit_truX1[ rep_i,i_loci  ] <- my_herit_lmm_gcta(name_pheno_X1, name_kinship_true)
    
    herit_truX2[ rep_i,i_loci  ] <- my_herit_lmm_gcta(name_pheno_X2, name_kinship_true)
    

    ###############
    ### CLEANUP ###
    ###############
    
    # these get deleted after every iteration
    #delete_files_plink(name_out)
    delete_files_phen(name_pheno_X1)
    delete_files_phen(name_pheno_X2)

    
    # herit output and log files (only two, they get overwritten for every GRM tested)
    delete_files_gcta(name_pheno_X1, herit = TRUE)
    delete_files_gcta(name_pheno_X2, herit = TRUE)
    
}
    i_loci = i_loci+1

}

###############
### CLEANUP ###
###############

# delete these shared files only after all replicates are done
delete_files_grm(name_kinship_true)


######################
### SAVE ESTIMATES ###
######################

# gather vectors into a dataframe
data <- cbind(
    herit_truX1[,1],
    herit_truX2[,1],
    herit_truX1[,2],
    herit_truX2[,2],
    herit_truX1[,3],
    herit_truX2[,3],
    herit_truN
)

# data <- cbind(
#     herit_truX1,
#     herit_truX2,
#     herit_truN
# )


herit_all1 <- as.data.frame(data)

names(herit_all1) <- c("true_genetic_1k_cst","true_genetic_1k","true_genetic_10k_cst","true_genetic_10k",
                       "true_genetic_100k_cst","true_genetic_100k","true_MVN")


write.csv(herit_all1,file = "herit_1k_10k_100k_v.3.csv",row.names = FALSE)

# setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/herit_1201/results_mvn_genetic_bias' )
# herit_all1<-read.csv("herit_100_1k_10k_100k_200K.csv")


###############
### boxplot ###
###############


# library(tidyverse)
# 
# herit_all1 %>% 
#     pivot_longer(names(herit_all1), names_to = "methods", values_to = "herit") ->herit_all2
# 
# herit_all2$methods <- factor(herit_all2$methods, levels =names(herit_all1))
# 
# ggplot(herit_all2, aes(x = methods, y = herit)) +
#     geom_boxplot() + labs(y = "heritability", x = "")+
#     geom_hline(yintercept=herit, linetype="dashed", color = "red")









