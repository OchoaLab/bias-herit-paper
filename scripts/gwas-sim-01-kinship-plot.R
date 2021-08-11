library(optparse) # for terminal options
library(genio)    # to read GRM files
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF

############
### ARGV ###
############

# define options
option_list = list(
    make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
                help = "number of individuals", metavar = "int"),
    make_option(c("-m", "--m_loci"), type = "integer", default = 10000, 
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
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model")
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
fes <- opt$fes

# output path for BED files and all results files
dir_out <- paste0(
    'sim-admix',
    '-n', n_ind,
    '-m', m_loci,
    '-k', k_subpops,
    '-f', fst,
    '-s', bias_coeff,
    '-mc', m_causal,
    '-h', herit,
    '-g', generations,
    if ( fes ) '-fes' else '-rc'
)

# load pre-existing data
setwd( '../data/' )
#setwd( 'D:/3.Duke/research/alex_ochoa/1.reverse_regression/coding/mycode/true-vs-biased-kinship-gwas/results' )
setwd( dir_out )

# they have these names
# all are 2x, scaled as GCTA wants them, halve here for plots
kinship_true <- read_grm( 'kinship_true' )$kinship / 2
kinship_popkin <- read_grm( 'kinship_popkin' )$kinship / 2
kinship_std_rom <- read_grm( 'kinship_std_rom' )$kinship / 2
kinship_std_rom_lim <- read_grm( 'kinship_std_rom_lim' )$kinship / 2
kinship_std_mor <- read_grm( 'kinship_std_mor' )$kinship / 2
kinship_gcta <- read_grm( 'kinship_gcta' )$kinship / 2
kinship_gcta_lim <- read_grm( 'kinship_gcta_lim' )$kinship / 2
kinship_wg <- read_grm( 'kinship_wg' )$kinship / 2
kinship_wg_lim <- read_grm( 'kinship_wg_lim' )$kinship / 2

# visualize all matrices for test
dims <- fig_scale( ratio = 3/4 )
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2]
)
plot_popkin(
    inbr_diag(
        list(
            kinship_true,
            kinship_popkin,
            NULL,
            kinship_std_rom_lim,
            kinship_std_rom,
            kinship_std_mor,
            kinship_gcta_lim,
            kinship_gcta,
            NULL,
            kinship_wg_lim,
            kinship_wg
        )
    ),
    titles = c(
        'Truth',
        'Popkin est.',
        'Standard ROM lim.',
        'Standard ROM est.',
        'Standard MOR est.',
        'GCTA lim.',
        'GCTA est.',
        'Weir-Goudet lim.',
        'Weir-Goudet est.'
    ),
    layout_rows = 4,
    mar = c(0, 2),
    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
