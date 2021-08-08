# simulates kinship matrices to test a hypothesis about the commutability of centering and EVD

# this function is for generating random genotypes from an admixed population
# defaults for 3 source subpopulations, meant to abstractly resemble Hispanics and African-Americans

# things that can't be changed:
# - Fst of intermediate subpopulations is a ramp proportional to 1:k
# - Admixture proportions are from 1D geography model

library(optparse) # for terminal options
library(bnpsd)    # simulate admixed population (structured population)
library(popkinsuppl) # for centering true kinship matrix (limit of biased standard kinship estimator)
library(popkin)   # plot_popkin
library(Matrix)   # rankMatrix
library(simfam)   # for simulating family structures

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-k", "--k_subpops"), type = "integer", default = 3, 
                help = "admixture intermediate subpopulations", metavar = "int"),
    make_option(c("-f", "--fst"), type = "double", default = 0.3, 
                help = "FST (fixation index)", metavar = "double"),
    make_option(c("--bias_coeff"), type = "double", default = 0.5, 
                help = "admixture bias coeff", metavar = "double"),
    make_option(c("-g", "--generations"), type = "integer", default = 1, 
                help = "number of generations, for realistic local kinship", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
G <- opt$generations

# hardcoded params
verbose <- TRUE # to print messages

############################
### POPULATION STRUCTURE ###
############################

# these are steps involving admixture and family structure abstractly, without genotypes at all (or other allele frequencies)

# define population structure
inbr_subpops <- 1 : k_subpops # FST values for 3 subpopulations
if (verbose)
    message('admix_prop_1d_linear')
obj <- admix_prop_1d_linear(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = fst
)
# in this case return value is a named list with three items:
admix_proportions <- obj$admix_proportions # admixture proportions
inbr_subpops <- obj$coanc_subpops # rescaled Fst vector for intermediate subpops

# calculate true kinship matrix
if (verbose)
    message('coanc_admix, coanc_to_kinship')
coancestry <- coanc_admix(admix_proportions, inbr_subpops) # the coancestry matrix
kinship <- coanc_to_kinship(coancestry) # kinship matrix

if ( G > 1 ) {
    # simulate realistic generations of families

    # simulate pedigree first
    if (verbose)
        message('sim_pedigree')
    data_simfam <- sim_pedigree( n_ind, G )
    fam <- data_simfam$fam # the pedigree itself
    ids <- data_simfam$ids # to filter generations later

    if (verbose)
        message('kinship_fam')
    # label founders in matrix
    rownames( kinship ) <- ids[[ 1 ]]
    colnames( kinship ) <- ids[[ 1 ]]
    # now calculate total kinship in final generation
    kinship <- kinship_last_gen( kinship, fam, ids )
}

# now test whether centering and low-dimensional approximation commute

# wrapper for EVD-based approximation
evd_approx <- function( kinship, r ) {
    # decompose it...
    evd <- eigen(kinship)
    # extract submatrices
    L <- evd$values
    U <- evd$vectors
    # indexes to keep
    indexes <- 1 : r
    # subset as required
    L <- L[ indexes ]
    U <- U[ , indexes ] # subset columns
    # reform kinship, the low-dimensional approximation
    kinship_r <- tcrossprod( U %*% diag(L), U )
    # return the kinship matrix and other data
    return(
        list(
            kinship = kinship_r,
            values = L,
            vectors = U
        )
    )
}

# first center, then EVD-approximate (this is what is done in practice, the limit of that anyway)
kinship_std_lim <- kinship_std_limit( kinship )
kinship_std_lim_r <- evd_approx( kinship_std_lim, k_subpops )$kinship

# now the other way, EVD first, then center that
kinship_r <- evd_approx( kinship, k_subpops )$kinship
kinship_r_std_lim <- kinship_std_limit( kinship_r )

## # visualize them all
## plot_popkin(
##     inbr_diag(
##         list(
##             kinship,
##             kinship_r,
##             kinship_std_lim,
##             kinship_std_lim_r,
##             kinship_r_std_lim
##         )
##     ),
##     titles = c(
##         'True kinship',
##         'True kinship, dim-r',
##         'Biased limit',
##         'Biased limit, dim-r',
##         'Dim-r, centered'
##     )
## )

## # plot the two of interest separately
## plot_popkin(
##     inbr_diag(
##         list(
##             kinship_std_lim_r,
##             kinship_r_std_lim
##         )
##     ),
##     titles = c(
##         'Biased limit, dim-r',
##         'Dim-r, centered'
##     )
## )

## # dotplots, for more precise comparisons
## # I should separate diagonals in this case (most likely source of disagreement), but meh
## plot(
##     inbr_diag(kinship), # transforming this but not the other agrees best here!
##     kinship_r,
##     xlab = 'True kinship',
##     ylab = 'True kinship, dim-r',
##     pch = '.'
## )
## abline( 0, 1, lty = 2, col = 'gray')

## # skip to the final two... the key comparison
## plot(
##     kinship_std_lim_r,
##     kinship_r_std_lim,
##     xlab = 'Biased limit, dim-r',
##     ylab = 'Dim-r, centered',
##     pch = '.'
## )
## abline( 0, 1, lty = 2, col = 'gray')

# a simple, numerical summary of agreement
message(
    'RMSD center/dim-r, vs dim-r/center: ',
    sqrt( mean( ( kinship_std_lim_r - kinship_r_std_lim )^2 ) )
)
# turn to vectors, then compute correlation coefficient
message(
    'Corr center/dim-r, vs dim-r/center: ',
    cor( as.numeric( kinship_std_lim_r ), as.numeric( kinship_r_std_lim ) )
)
# sum of kinship, should be zero if the intercept is in the nullspace
message(
    'Sum of kinship, center/dim-r: ',
    signif( sum( kinship_std_lim_r ), 3 )
)
message(
    'Sum of kinship, dim-r/center: ',
    signif( sum( kinship_r_std_lim ), 3 )
)
# direct measurement of matrix rank, including verification of overlapping rowspaces
message(
    'Rank center/dim-r: ',
    rankMatrix( kinship_std_lim_r )
)
message(
    'Rank dim-r/center: ',
    rankMatrix( kinship_r_std_lim )
)
message(
    'Rank center/dim-r + dim-r/center: ',
    rankMatrix( cbind( kinship_std_lim_r, kinship_r_std_lim ) )
)
message(
    'Rank center/dim-r + dim-r: ',
    rankMatrix( cbind( kinship_std_lim_r, kinship_r ) )
)
message(
    'Rank dim-r: ',
    rankMatrix( kinship_r )
)
message(
    'Rank dim-r/center: ',
    rankMatrix( kinship_r_std_lim )
)
# the actual rowspaces of interest include the intercept!
message(
    'Rank intercept + dim-r: ',
    rankMatrix( cbind( kinship_r, 1 ) )
)
message(
    'Rank intercept + dim-r/center: ',
    rankMatrix( cbind( kinship_r_std_lim, 1 ) )
)
message(
    'Rank intercept + center/dim-r: ',
    rankMatrix( cbind( kinship_std_lim_r, 1 ) )
)
message(
    'Rank intercept + dim-r/center + dim-r: ',
    rankMatrix( cbind( kinship_r_std_lim, kinship_r, 1 ) )
)
message(
    'Rank intercept + center/dim-r + dim-r: ',
    rankMatrix( cbind( kinship_std_lim_r, kinship_r, 1 ) )
)
