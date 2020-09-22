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

# switch to main scripts directory
setwd( '../../scripts/' )


# standard code for a complex trait and an admixed population
source('sim_geno_trait_k3.R')

# load new functions from external scripts
source('kinship_to_evd.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')

# get back to this subdirectory
setwd( '../bias/scripts/' )


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

############
### GCTA ###
############

message( "GCTA (default)" )
# estimate kinship with GCTA's default methodd
gas_lmm_gcta_kin(gcta_bin, name_out)
# GWAS using that default kinship matrix
obj_gcta_default <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads)
# obj_gcta_default$pvals
# obj_gcta_default$beta_hat
# cleanup
delete_files_gcta(name_out) # GAS table
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
obj_gcta_true <- gas_lmm_gcta(
    gcta_bin,
    name_out,
    name_grm = name_kinship_true,
    m_loci = m_loci,
    threads = threads
)
## obj_gcta_true$pvals
## obj_gcta_true$beta_hat
# cleanup
delete_files_gcta(name_out) # GAS table
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
obj_gcta_popkin <- gas_lmm_gcta(
    gcta_bin,
    name_out,
    name_grm = name_kinship_popkin,
    m_loci = m_loci,
    threads = threads
)
## obj_gcta_popkin$pvals
## obj_gcta_popkin$beta_hat
# cleanup
delete_files_gcta(name_out) # GAS table
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
obj_gcta_std_lim <- gas_lmm_gcta(
    gcta_bin,
    name_out,
    name_grm = name_kinship_std_lim,
    m_loci = m_loci,
    threads = threads
)
## obj_gcta_std_lim$pvals
## obj_gcta_std_lim$beta_hat
# cleanup
delete_files_gcta(name_out) # GAS table
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
## obj_gcta_std$pvals
## obj_gcta_std$beta_hat
# cleanup
delete_files_gcta(name_out) # GAS table
delete_files_grm(name_kinship_std) # GRM

# final cleanup: do after all GCTA runs are finished
# delete the temporary BED files now that we are done
delete_files_plink(name_out)
delete_files_phen(name_out)

#######

gcta_pval<-data.frame(obj_gcta_default$pvals,obj_gcta_true$pvals,obj_gcta_popkin$pvals,obj_gcta_std_lim$pvals,obj_gcta_std$pvals)

colnames(gcta_pval)<-c("gcta_default","gcta_true","gcta_popkin","gcta_std_lim","gcta_std")

cor(gcta_pval)

#write.csv(gcta_pval)

gcta_beta_hat<-data.frame(obj_gcta_default$beta_hat,obj_gcta_true$beta_hat,obj_gcta_popkin$beta_hat,obj_gcta_std_lim$beta_hat,obj_gcta_std$beta_hat)

colnames(gcta_beta_hat)<-c("gcta_default","gcta_true","gcta_popkin","gcta_std_lim","gcta_std")

cor(gcta_beta_hat)



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

######

pca_pval<-data.frame(obj_pca_true$pvals,obj_pca_popkin$pvals,obj_pca_std_lim$pvals,obj_pca_std$pvals)

colnames(pca_pval)<-c("pca_true","pca_popkin","pca_std_lim","pca_std")

cor(pca_pval)


pca_beta_hat<-data.frame(obj_pca_true$beta_hat,obj_pca_popkin$beta_hat,obj_pca_std_lim$beta_hat,obj_pca_std$beta_hat)

colnames(pca_beta_hat)<-c("pca_true","pca_popkin","pca_std_lim","pca_std")

cor(pca_beta_hat)

###############
### Results ###
###############


###### mixed
# pvaule

all_pval<-data.frame(gcta_pval,pca_pval)

cor(all_pval)

write.table(all_pval,"all_pval.txt",row.names = FALSE  )

#beta_hat

all_beta_hat<-data.frame(gcta_beta_hat,pca_beta_hat)

cor(all_beta_hat)

write.table(all_beta_hat,"all_beta_hat.txt",row.names = FALSE  )


###############################
### Heatmap & Clustering ----
###############################

library(dendextend)
# order for rows
Rowv  <-  all_pval %>%  t %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_lwd", 1.2) %>%
    ladderize

# Order for columns: We must transpose the data
Colv  <- all_pval %>%  t %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_lwd", 1.2) %>%
    ladderize


#heatmap(as.matrix(cor(all_pval)), Rowv = Rowv, Colv = Colv,scale = "none")


library(gplots)
library(RColorBrewer)
heatmap.2(as.matrix(cor(all_pval)), scale = "none",
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          breaks = seq(-1,1,0.02),
          Rowv = Rowv, Colv = Colv,trace = "none",cexRow =1,cexCol=0.9,
          offsetRow=0,offsetCol=0,density.info = "none",
          key = TRUE,key.title=NA)

heatmap.2(as.matrix(cor(all_beta_hat)), scale = "none",
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          breaks = seq(-1,1,0.02),
          Rowv = Rowv, Colv = Colv,trace = "none",cexRow =1,cexCol=0.9,
          offsetRow=0,offsetCol=0,density.info = "none",
          key = TRUE,key.title=NA)

# library(corrplot)
# corrplot(cor(all_pval), type="upper", order="hclust")




