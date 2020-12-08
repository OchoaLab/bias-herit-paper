# compare details of GWAS with true or biased kinship limits, and the two kinds of (noisy) estimates

library(optparse)    # for terminal options
library(readr)       # to write kinship matrix
library(tibble)      # to store data
library(genio)       # to write BED files for external software
library(popkin)      # to estimate kinship in RGLS
library(popkinsuppl) # for PCA's kinship estimator

# switch to main scripts directory
dir_orig <- getwd() # remember where we are now
#setwd( '../../scripts/' )
setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/gas-rgls/scripts' )

source('sim_geno_trait_k3.R')
source('gas_lmm_gcta.R')
source('kinship_gcta_limit.R')
#source('herit_lmmlite.R')
#source('paths.R')
setwd( dir_orig ) # go back to where we were

# now move to data location
setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/herit_1201/results_kinship_bias' )


# change this path to GCTA to whatever it is in your computer
#gcta_bin <- '/home/viiia/bin/gcta_1.93.2beta/gcta64'
gcta_bin <- 'D:/3.Duke/research/alex_ochoa/software/gcta_1.93.2beta_win/bin/gcta64'

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-r", "--rep"), type = "integer", default = 10, 
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
generations <- opt$generations
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
    '-g', generations,
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
    generations = generations,
    #    iterations = 100,
    verbose = TRUE
)
admix_proportions <- obj$admix_proportions
inbr_subpops <- obj$inbr_subpops
parents <- obj$parents
kinship <- obj$kinship

# There are 6 kinship matrices to consider
# 1- `kinship`: true kinship matrix of simulation
# 2- `kinship_popkin`: unbiased (but noisy) estimate from genotypes


# 1) true kinship (as given by simulation)


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

loci_test <- c(10000,100000,200000)

# store the heritability estimates in these vectors

# then using trait from MVN model
herit_truN <- vector( 'numeric', rep )


herit_popN <- matrix( NA, ncol = length(loci_test), nrow = rep )


#################################################
### Estimate heritability based on true kinship##
#################################################

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
### Estimate heritability based on popkin ########
##################################################

i_loci <-1

for ( m_loci in loci_test  ) {

    for ( rep_i in 1 : rep ) {
    
    ################
    ### SIM GENO ###
    ################
    
    # simulate new genotypes for each replicate
    obj <- sim_geno(
        admix_proportions = admix_proportions,
        inbr_subpops = inbr_subpops,
        generations = generations,
        parents = parents,
        m_loci = m_loci,
        verbose = TRUE
    )
    X <- obj$X
    
    # and write MVN trait too
    fam$pheno <- traits_mvn[ rep_i, ]
    name_pheno_N <- paste0(name_out, '.pheno_N')
    write_phen(name_pheno_N, fam)
    
    ###############
    ### KINSHIP ###
    ###############
    
    # estimates that get redone at every replicate
    
    # 2) popkin estimate
    # need labels first
    labs <- ceiling( ( 1 : n_ind ) / n_ind * 10 )
    # actual popkin estimate
    message( 'popkin' )
    kinship_popkin <- popkin(X, labs)
    

    # SAVE as GRM files
    name_kinship_popkin <- paste0(name_out, '.kinship_popkin')
    my_write_grm( name_kinship_popkin, kinship_popkin )
    

    
    ############
    ### GCTA ###
    ############
    
    # get estimates
    # another wrapper to avoid mistakes in these loops
    my_herit_lmm_gcta <- function ( name_pheno, name_kinship ) {
        herit_lmm_gcta(gcta_bin, name_pheno, name_grm = name_kinship, threads = threads)$herit
    }
    

    # then using trait from MVN model
    
    message('herit_popN...')
    herit_popN[ rep_i, i_loci ] <- my_herit_lmm_gcta(name_pheno_N, name_kinship_popkin)
    

    ###############
    ### CLEANUP ###
    ###############
    
    # these get deleted after every iteration
    delete_files_grm(name_kinship_popkin)

}
    i_loci = i_loci+1
    

}

# cleanup
delete_files_log(name_kinship_gcta)



# gather vectors into a dataframe
data <- cbind(
    herit_popN,
    herit_truN
)

herit_all1 <- as.data.frame(data)
names(herit_all1) <- c("popkin_10k","popkin_100k","popkin_200k",
                       "true_MVN")

#write.csv(herit_all1,file = "herit_200k_1207.csv",row.names = FALSE)


#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/heritability_estimation' )


#herit_all<-read.csv("herit_lmmlite_n=2000_50times-raw.csv")



###############
### boxplot ###
###############



library(tidyverse)

herit_all1 %>% 
    pivot_longer(names(herit_all1), names_to = "methods", values_to = "herit") ->herit_all2

herit_all2$methods <- factor(herit_all2$methods, levels =names(herit_all1))

ggplot(herit_all2, aes(x = methods, y = herit)) +
    geom_boxplot() + labs(y = "heritability", x = "")+
    geom_hline(yintercept=herit, linetype="dashed", color = "red")



