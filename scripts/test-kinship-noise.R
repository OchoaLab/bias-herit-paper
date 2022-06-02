# confirms heritability biases due to noisy estimation of otherwise unbiased kinship, reduced by increasing the number of loci

library(optparse)    # for terminal options
library(readr)       # to write kinship matrix
library(genio)       # to write BED files for external software
library(popkin)      # to estimate kinship in RGLS
library(simgenphen)  # for simulating genotypes and phenotypes
library(genbin)      # to estimate heritability with binary GCTA
library(simtrait)    # for sim_trait_mvn
library(ochoalabtools) # for fig

# now move to data location
setwd( '../data' )

############
### ARGV ###
############

# hardcoded numbers of loci to simulate
m_loci_test <- c(10000,100000,200000)

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
G <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
threads <- opt$threads

# output path for BED files
name_out <- paste0(
    'herit-kinship-noise',
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

#############
### TESTS ###
#############

# store the heritability estimates in these vectors
herit_truN <- vector( 'numeric', rep )
herit_popN <- matrix( NA, ncol = length(m_loci_test), nrow = rep )

for ( rep_i in 1 : rep ) {
    # draw population structure (differs across replicates if G > 1)
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

    # SAVE true kinship as GRM (shared across replicates)
    # matrix of sample sizes is complete matrix in all these cases (no missingness)
    M <- matrix( m_loci, nrow = n_ind, ncol = n_ind )
    # also need a dummy fam file for writing the GRM, same used for phenotypes (so whatever we do they match!)
    fam <- make_fam( n = n_ind )
    # write the kinship matrix to a file
    name_kinship_true <- paste0(name_out, '.kinship_true')
    write_grm(
        name_kinship_true,
        kinship = 2 * kinship,
        M = M,
        fam = fam
    )

    message( 'sim_trait_mvn...' )
    # draw single MVN trait
    # (although drawing multiple traits is more efficient for fixed kinship, here kinship may not be fixed if G > 1)
    fam$pheno <- drop( sim_trait_mvn(
        rep = 1,
        kinship = kinship,
        herit = herit
    ) )
    # write trait to file (MVN is only case tested here)
    name_pheno_N <- paste0( name_out, '.pheno_N' )
    write_phen( name_pheno_N, fam )
    
    # get herit estimate with true kinship
    message('herit_truN...')
    herit_truN[ rep_i ] <- gcta_reml( name_pheno_N, name_grm = name_kinship_true, threads = threads )$herit
    # cleanup
    delete_files_gcta_hsq( name_pheno_N )
    delete_files_log( name_pheno_N )
    delete_files_grm( name_kinship_true )

    for ( i_loci in 1 : length( m_loci_test ) ) {
        # simulate new genotypes for each replicate
        obj <- sim_geno(
            admix_proportions_1 = admix_proportions_1,
            inbr_subpops = inbr_subpops,
            fam = famG,
            ids = ids,
            m_loci = m_loci_test[ i_loci ]
        )
        X <- obj$X

        # estimates that get redone at every replicate
        message( 'popkin' )
        kinship_popkin <- popkin( X )
        
        # SAVE popkin estimate as GRM (changes for every rep and m_loci)
        name_kinship_popkin <- paste0(name_out, '.kinship_popkin')
        write_grm(
            name_kinship_popkin,
            kinship = 2 * kinship_popkin,
            M = M,
            fam = fam
        )

        # get estimates
        # then using trait from MVN model
        message('herit_popN...')
        herit_popN[ rep_i, i_loci ] <- gcta_reml( name_pheno_N, name_grm = name_kinship_popkin, threads = threads )$herit
        # cleanup
        delete_files_gcta_hsq( name_pheno_N )
        delete_files_log( name_pheno_N )
        delete_files_grm( name_kinship_popkin )
    }
    
    # final rep cleanup
    delete_files_phen( name_pheno_N )
}


##################
### SAVE TABLE ###
##################

# gather vectors into a dataframe
data <- cbind(
    herit_popN,
    herit_truN
)

data <- as.data.frame(data)
names(data) <- c( "popkin_10k", "popkin_100k", "popkin_200k", "true" )

write_tsv( data, file = paste0( name_out, ".txt" ) )

###############
### boxplot ###
###############

library(tidyverse)

data %>% 
    pivot_longer(names(data), names_to = "methods", values_to = "herit") -> data2

data2$methods <- factor(data2$methods, levels =names(data))

fig_start( name_out, width = 6 )

ggplot(data2, aes(x = methods, y = herit)) +
    geom_boxplot() + labs(y = "heritability", x = "")+
    geom_hline(yintercept=herit, linetype="dashed", color = "red")

fig_end()

