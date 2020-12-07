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
setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/herit_1201/results' )


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



###############
### KINSHIP ###
###############

# There are 2 kinship matrices to consider
# - `kinship`: true kinship matrix of simulation
# - `kinship_popkin`: unbiased (but noisy) estimate from genotypes

# 1) true kinship (as given by simulation)

# 2) popkin estimate
# need labels first
labs <- ceiling( ( 1 : n_ind ) / n_ind * 10 )
# actual popkin estimate
kinship_popkin <- popkin(X, labs)


#################
### SIM TRAIT ###
#################

rep <-10

message( 'sim_trait_mvn...' )

# draw all traits at once (most efficient)
traits <- sim_trait_mvn(
    rep = rep,
    kinship = kinship,
    herit = herit
)

#######


# initialize result matrix
herit_all <- matrix( NA, ncol = 2, nrow = rep )

for ( i in 1 : rep ) {
    trait <- traits[ i, ]
    

    # write BED version for external code
    plink_data <- write_plink(name_out, X, pheno = traits[i,], verbose = FALSE)
    # write phenotype file
    write_phen(name_out, plink_data$fam, verbose = FALSE)

    # matrix of sample sizes is complete matrix in all these cases (no missingness)
    M <- matrix( m_loci, nrow = n_ind, ncol = n_ind )
    
    message( "GCTA (true kinship)" )
    # write the kinship matrix to a file
    name_kinship_true <- paste0(name_out, '.kinship_true')
    write_grm(
        name_kinship_true,
        kinship = 2 * kinship,
        M = M,
        fam = plink_data$fam,
        verbose = FALSE
    )

    obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_true, threads = threads)
    herit_all[i,1]<-obj1$herit
    rm(obj1)

    # cleanup
    #delete_files_gcta(name_out) # GAS table
    delete_files_grm(name_kinship_true) # GRM


    message( "GCTA (popkin kinship)" )
    # write the kinship matrix to a file
    name_kinship_popkin <- paste0(name_out, '.kinship_popkin')
    write_grm(
        name_kinship_popkin,
        kinship = 2 * kinship_popkin,
        M = M,
        fam = plink_data$fam,
        verbose = FALSE
    )

    obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_popkin, threads = threads)
    herit_all[i,2]<-obj1$herit
    rm(obj1)


    # cleanup
    #delete_files_gcta(name_out) # GAS table
    delete_files_grm(name_kinship_popkin) # GRM

    #final cleanup: do after all GCTA runs are finished
    #delete the temporary BED files now that we are done
    delete_files_plink(name_out)
    delete_files_phen(name_out)

}

# rename results matrix 

#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/heritability_estimation' )
herit_all_200k <- as.data.frame(herit_all)
names(herit_all_200k) <- c("true kin_200k","popkin_200k")


herit_all_100k <- as.data.frame(herit_all)
names(herit_all_100k) <- c("true kin_100k","popkin_100k")


herit_all_10k <- as.data.frame(herit_all)
names(herit_all_10k) <- c("true kin_10k","popkin_10k")

herit_all1<-cbind(herit_all_200k,herit_all_100k,herit_all_10k)

#write.csv(herit_all1,file = "herit_200k.csv",row.names = FALSE)


#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/heritability_estimation' )


#herit_all<-read.csv("herit_lmmlite_n=2000_50times-raw.csv")

#names(herit_all)<-c("default","true kin","popkin","std lim kin","std kin","gcta lim kin")


## boxplot

library(tidyverse)

herit_all1 %>% 
    pivot_longer(c("true kin_200k","popkin_200k","true kin_100k","popkin_100k",
                   "true kin_10k","popkin_10k"), names_to = "methods", values_to = "herit") ->herit_all2


herit_all2$methods <- factor(herit_all2$methods, levels =
                                 c("true kin_200k","popkin_200k","true kin_100k","popkin_100k",
                                   "true kin_10k","popkin_10k"))

#factor(herit_all2$methods)

herit_all2$methods<-as.factor(herit_all2$methods)

ggplot(herit_all2, aes(x = methods, y = herit)) +
    geom_boxplot() + labs(y = "heritability", x = "")+
    geom_hline(yintercept=herit, linetype="dashed", color = "red")



save.image(file = "herit_200k.RData")


