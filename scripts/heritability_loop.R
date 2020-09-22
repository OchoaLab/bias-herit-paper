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

# get back to this subdirectory
setwd( '../env/scripts/' )


# change this path to GCTA to whatever it is in your computer
gcta_bin <- '/home/viiia/bin/gcta_1.93.2beta/gcta64'



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
  make_option(c("-t", "--threads"), type = "integer", default = 1, 
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
  '-g', generations
)

# initialize result matrix

n<-10

herit_all<-matrix(NA,ncol = 6, nrow = n )



for (i in 1:n) {


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
# write phenotype file
write_phen(name_out, plink_data$fam, verbose = FALSE)


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


# 5) gcta limit estimate
kinship_gcta_lim <- kinship_gcta_limit( kinship )


############
### GCTA ###
############

message( "GCTA (default)" )


# estimate kinship with GCTA's default methodd
gas_lmm_gcta_kin(gcta_bin, name_out)

# GWAS using that default kinship matrix
#obj_gcta_default <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads)

obj1 <- herit_lmm_gcta(gcta_bin, name_out, threads = threads)
herit_all[i,1]<-obj1$herit
rm(obj1)

# cleanup
#delete_files_gcta(name_out) # GAS table
delete_files_grm(name_out) # GRM default method



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
# GWAS
# obj_gcta_true <- gas_lmm_gcta(
#   gcta_bin,
#   name_out,
#   name_grm = name_kinship_true,
#   m_loci = m_loci,
#   threads = threads
# )

obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_true, threads = threads)
herit_all[i,2]<-obj1$herit
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
# GWAS
# obj_gcta_popkin <- gas_lmm_gcta(
#   gcta_bin,
#   name_out,
#   name_grm = name_kinship_popkin,
#   m_loci = m_loci,
#   threads = threads
# )


obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_popkin, threads = threads)
herit_all[i,3]<-obj1$herit
rm(obj1)


# cleanup
#delete_files_gcta(name_out) # GAS table
delete_files_grm(name_kinship_popkin) # GRM


# GCTA with limit of biased "standard" kinship estimate
message( "GCTA (std lim kinship)" )
# write the kinship matrix to a file
name_kinship_std_lim <- paste0(name_out, '.kinship_std_lim')
write_grm(
  name_kinship_std_lim,
  kinship = 2 * kinship_std_lim,
  M = M,
  fam = plink_data$fam,
  verbose = FALSE
)
# GWAS
# obj_gcta_std_lim <- gas_lmm_gcta(
#   gcta_bin,
#   name_out,
#   name_grm = name_kinship_std_lim,
#   m_loci = m_loci,
#   threads = threads
# )

obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_std_lim, threads = threads)
herit_all[i,4]<-obj1$herit
rm(obj1)



# cleanup
#delete_files_gcta(name_out) # GAS table
delete_files_grm(name_kinship_std_lim) # GRM





# GCTA with biased "standard" kinship estimate
message( "GCTA (std kinship)" )
# write the kinship matrix to a file
name_kinship_std <- paste0(name_out, '.kinship_std')
write_grm(
  name_kinship_std,
  kinship = 2 * kinship_std,
  M = M,
  fam = plink_data$fam,
  verbose = FALSE
)
# GWAS
obj_gcta_std <- gas_lmm_gcta(
  gcta_bin,
  name_out,
  name_grm = name_kinship_std,
  m_loci = m_loci,
  threads = threads
)


obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_std, threads = threads)
herit_all[i,5]<-obj1$herit
rm(obj1)



# cleanup
#delete_files_gcta(name_out) # GAS table
delete_files_grm(name_kinship_std) # GRM


message( "GCTA (limit)" )

name_kinship_gcta_lim <- paste0(name_out, '.kinship_gcta_lim')
write_grm(
  name_kinship_gcta_lim,
  kinship = 2 * kinship_gcta_lim,
  M = M,
  fam = plink_data$fam,
  verbose = FALSE
)

obj1 <- herit_lmm_gcta(gcta_bin, name_out,name_grm = name_kinship_gcta_lim, threads = threads)
herit_all[i,6]<-obj1$herit
rm(obj1)

delete_files_grm(name_kinship_gcta_lim) # GRM





# final cleanup: do after all GCTA runs are finished
# delete the temporary BED files now that we are done
delete_files_plink(name_out)
delete_files_phen(name_out)

}

## rename results matrix 

setwd( '../env/scripts/' )

herit_all<-as.data.frame(herit_all)

names(herit_all)<-c("default","true kin","popkin","std lim kin","std kin","gcta lim kin")

herit_all

#write.csv(herit_all,file = "herit_all_6_n_ind=2000-raw.csv",row.names = FALSE)


herit_all<-read.csv("herit_all_6_n_ind=2000-raw.csv")

names(herit_all)<-c("default","true kin","popkin","std lim kin","std kin","gcta lim kin")


## boxplot/violin plot

library(tidyverse)

herit_all %>% 
  pivot_longer(c("default","true kin","popkin","std lim kin","std kin","gcta lim kin"), names_to = "methods", values_to = "herit") ->herit_all2


herit_all2$methods <- factor(herit_all2$methods, levels =
                               c("default","true kin","popkin","std lim kin","std kin","gcta lim kin"))

herit_all2$methods<-as.factor(herit_all2$methods)

ggplot(herit_all2, aes(x = methods, y = herit)) +
  geom_boxplot() + labs(y = "heritability", x = "")+
  geom_hline(yintercept=0.8, linetype="dashed", color = "red")







