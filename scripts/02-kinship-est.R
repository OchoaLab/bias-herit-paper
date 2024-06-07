# simulate genotypes and phenotypes in reps

library(optparse)    # for terminal options
library(genio)       # for write_grm
library(BEDMatrix)   # for input to popkin
library(popkin)      # for popkin
library(popkinsuppl) # for kinship_std estimator
library(genbin)      # for running external binaries

# constants
# output data names (genotypes and trait)
name_out <- 'data'

# define options
option_list = list(
    make_option("--name", type = "character", default = NA, 
                help = "Base name for genotype and phenotype simulation", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "replicate number", metavar = "int"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (affects GCTA only)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$name
rep <- opt$rep
threads <- opt$threads

# now move to data location
setwd( '../data/' )
# move into project
setwd( dir_out )
# dir must already exist since we need genotypes present!
setwd( paste0( 'rep-', rep ) )

# helper functions

# this version overwrites new file if it existed
softlink <- function( name1, name2, ext ) {
    file1 <- paste0( name1, '.', ext )
    file2 <- paste0( name2, '.', ext )
    system3( 'ln', c( '-f', '-s', file1, file2 ) )
}

# automated version for most cases
write_grm_softlinks <- function(name, kinship) {
    # write main GRM files
    write_grm( name, 2 * kinship )
    # create a softlink from true "M" file, so those are the same for all files (limit variation due to that)
    # same with ID table (less big but whatever)
    # softlink saves space!
    softlink( 'std_mor', name, 'grm.id' ) # these files exist but are bad, overwrite!
    softlink( 'std_mor', name, 'grm.N.bin' )
}

softlink_limits <- function( name ) {
    # link precomputed limit kinship from base of project (shared across reps)
    softlink( paste0( '../../../kinship/', name ), name, 'grm.bin' )
    # these are same as local std_mor, as is for all other cases
    softlink( 'std_mor', name, 'grm.id' )
    softlink( 'std_mor', name, 'grm.N.bin' )
}

# put all kinship matrices in a subdirectory, which hasn't been created!
dir.create( 'kinship' )

# standard estimate, calculated with GCTA
# (only case not calculated in R, need gcta64 binary)
# do first as other estimators link aux files to these
gcta_grm( name_out, name_out = 'kinship/std_mor', threads = threads )
delete_files_log( 'kinship/std_mor' ) # unneeded log file

# load genotypes for kinship estimates
X <- BEDMatrix( name_out, simple_names = TRUE )

setwd( 'kinship' )

# for additional ease, add links of limits shared across reps
softlink_limits( 'true' )
softlink_limits( 'std_rom_lim' )


# biased "standard" kinship estimate, ROM only (because MOR == GCTA)
write_grm_softlinks( 'std_rom', kinship_std( X ) )

# popkin estimates, without labels
kinship_popkin_rom <- popkin( X )
write_grm_softlinks( 'popkin_rom', kinship_popkin_rom )
kinship_popkin_mor <- popkin( X, mean_of_ratios = TRUE )
write_grm_softlinks( 'popkin_mor', kinship_popkin_mor )


