# this script generates a few replicates of the same simulation with slight differences, to showcase the biases Zhuoran has identified.
# this follow up actually creates the plot

library(optparse)    # for terminal options
library(readr)       # to write data
library(ochoalabtools) # for plotting
library(tibble) # for a mapping table

# move to data location
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
    make_option("--grant", action = "store_true", default = FALSE, 
                help = "filter for grant (MIRA)")
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
grant <- opt$grant

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

# load data
data <- read_tsv(
    paste0( name_out, '.txt' )
)

# a map from internal codes to nice names
# first most straightforward map
codes <- c('tru', 'pop', 'slm', 'std', 'glm', 'gct')
nices <- c('Truth', 'Popkin', 'Standard Lim.', 'Standard', 'GCTA Lim.', 'GCTA')
# then doubled map for actual data (every estimator is tested on two kinds of traits)
methods <- tibble(
    code = c(
        paste0( codes, 'X' ),
        paste0( codes, 'N' )
    ),
    name = c(nices, nices),
    oracle = rep.int( c(1, 0), 6 ) # defined and doubled here
)

# difference for grant is mostly filtering data
if ( grant ) {
    # change output path
    name_out <- paste0( name_out, '-grant' )

    # remove data that would be confusing to show in grant (std_lim and gcta_lim, names are weirder though, meh)
    # list of things to keep
    names_keep <- c('tru', 'pop', 'std', 'gct')
    # add X and N suffixes
    names_keep <- c(
        paste0( names_keep, 'X' ),
        paste0( names_keep, 'N' )
    )
    # apply filter
    data <- data[ , names(data) %in% names_keep ]
}


# explicitly map data to labels in table, so these always agree even if there is reordering or filtering in the data
# filter and reorder methods
methods <- methods[ match( names(data), methods$code ), ]

col_oracle <- 'gray50'
col_estimated <- 'darkblue'

bot_label_line <- if (grant) 4 else 6

# begin plot
fig_start(
    name_out,
    width = if (grant) 3 else 4.5,
    height = if (grant) 3.5 else 4,
    mar_t = 2,
    mar_b = bot_label_line + 1
)
boxplot(
    data,
    names = methods$name,
    xlab = '',
    ylab = 'Heritability estimate',
    las = 3,
    border = if (grant) 'black' else ifelse( methods$oracle, col_oracle, col_estimated ),
    col = if (grant) 'lightgray' else 'lightgray'
)
if (!grant) {
    legend(
        if ( generations == 1 ) 'bottomright' else 'topright',
        title = 'Kinship type',
        title.col = 'black', # default uses col_oracle for some reason :(
        legend = c('Oracle', 'Estimated'),
        text.col = c(col_oracle, col_estimated),
        col = NA,
#        bty = 'n',
        cex = 0.7
        )
}
abline( h = herit, lty = 2, col = 'red' )
text(
    x = if ( generations == 1 ) ncol(data) + 0.5 else 0.5,
    y = herit,
    labels = "True\nHeritability",
    col = 'red',
    adj = if ( generations == 1 ) 1 else 0
)
mtext(
    'Kinship estimator',
    side = 1,
    line = bot_label_line
)
mtext(
    'Trait simulation type',
    side = 3,
    line = 1
)
mtext(
    'Genetic',
    adj = if (grant) 0.15 else 0.2
)
mtext(
    'MVN',
    adj = if (grant) 0.8 else 0.77
)
abline( v = ( ncol(data) + 1 ) / 2 )
fig_end()
# and trait type
