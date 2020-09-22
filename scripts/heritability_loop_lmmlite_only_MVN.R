# compare details of GWAS with true or biased kinship limits, and the two kinds of (noisy) estimates

library(optparse)    # for terminal options
library(readr)       # to write kinship matrix
library(tibble)      # to store data
library(genio)       # to write BED files for external software
library(popkin)      # to estimate kinship in RGLS
library(popkinsuppl) # for PCA's kinship estimator

# install.packages("devtools") # if needed
# library(devtools)
# install_github("OchoaLab/genio", build_opts = c())
# 

# switch to main scripts directory
setwd( '../../scripts/' )

# standard code for a complex trait and an admixed population
source('sim_geno_trait_k3.R')

# load new functions from external scripts
source('kinship_to_evd.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')
source('kinship_gcta_limit.R')

source('herit_lmmlite.R')

# get back to this subdirectory
setwd( '../env/scripts/' )

# change this path to GCTA to whatever it is in your computer
#gcta_bin <- '/home/viiia/bin/gcta_1.93.2beta/gcta64'

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
    make_option(c("-t", "--threads"), type = "integer", default = 1, 
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
    'gas',
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

################
### SIM GENO ###
################

# simulate kinship and genotypes only once, shared across replicates
# only trait will vary in each replicate

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

# # write BED version for external code
# plink_data <- write_plink(name_out, X, pheno = trait, verbose = FALSE)
# # write phenotype file
# write_phen(name_out, plink_data$fam, verbose = FALSE)

###############
### KINSHIP ###
###############

# There are 4 kinship matrices to consider
# - `kinship`: true kinship matrix of simulation
# - `kinship_std_lim`: limit of biased "standard" estimator
# - `kinship_std`: biased (and noisy) "standard" estimate from genotypes (same as GEMMAs)
# - `kinship_popkin`: unbiased (but noisy) estimate from genotypes

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

# # 5) gcta limit estimate
# kinship_gcta_lim <- kinship_gcta_limit( kinship )

#########################
### LMMLITE PREROTATE ###
#########################

message( 'herit_lmmlite_prerotate...' )

# saves time factoring an eigendecomposition for LMMLITE
prerot_kinship <- herit_lmmlite_prerotate( kinship )
prerot_kinship_std <- herit_lmmlite_prerotate( kinship_std )
prerot_kinship_std_lim <- herit_lmmlite_prerotate( kinship_std_lim )
prerot_kinship_popkin <- herit_lmmlite_prerotate( kinship_popkin )

#################
### SIM TRAIT ###
#################

message( 'sim_trait_mvn...' )

# draw all traits at once (most efficient)
traits <- sim_trait_mvn(
    rep = rep,
    kinship = kinship,
    herit = herit
)

###############
### lmmlite ###
###############

# initialize result matrix
herit_all <- matrix( NA, ncol = 4, nrow = rep )

for ( i in 1 : rep ) {
    trait <- traits[ i, ]
    
    message( "lmmlite (true kinship)" )
    herit_all[i,1] <- herit_lmmlite(trait, prerot_kinship, prerotated = TRUE)

    message( "lmmlite (popkin kinship)" )
    herit_all[i,2] <- herit_lmmlite(trait, prerot_kinship_popkin, prerotated = TRUE)

    #lmmlite with limit of biased "standard" kinship estimate
    message( "lmmlite (std lim kinship)" )
    herit_all[i,3] <- herit_lmmlite(trait, prerot_kinship_std_lim, prerotated = TRUE)

    # lmmlite with biased "standard" kinship estimate
    message( "lmmlite (std kinship)" )
    herit_all[i,4] <- herit_lmmlite(trait, prerot_kinship_std, prerotated = TRUE)
}

# rename results matrix 

#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/heritability_estimation' )
herit_all <- as.data.frame(herit_all)
names(herit_all) <- c("lmmlite_true","lmmlite_popkin","lmmlite_std_lim_kin","lmmlite_std_kin")

herit_all

#write.csv(herit_all,file = "herit_lmmlite_n=2000_50times-raw.csv",row.names = FALSE)


#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/heritability_estimation' )


#herit_all<-read.csv("herit_lmmlite_n=2000_50times-raw.csv")

#names(herit_all)<-c("default","true kin","popkin","std lim kin","std kin","gcta lim kin")


## boxplot/violin plot

library(tidyverse)

herit_all %>% 
    pivot_longer(c("lmmlite_true","lmmlite_popkin","lmmlite_std_lim_kin","lmmlite_std_kin"), names_to = "methods", values_to = "herit") ->herit_all2


herit_all2$methods <- factor(herit_all2$methods, levels =
                                                     c("lmmlite_true","lmmlite_popkin","lmmlite_std_lim_kin","lmmlite_std_kin"))

#factor(herit_all2$methods)

herit_all2$methods<-as.factor(herit_all2$methods)

ggplot(herit_all2, aes(x = methods, y = herit)) +
    geom_boxplot() + labs(y = "heritability", x = "")+
    geom_hline(yintercept=herit, linetype="dashed", color = "red")






