
# Install Libraries -------------------------------------------


## Bioconductor Packages ---------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

for (pkg in c(
              # General Packages
                "knitr",
                "tinytex",
                "devtools",
                
              # Microbiome Packages
                "dada2", 
                "phyloseq", 
                "ALDEx2",
                "Maaslin2",
                "ANCOMBC",
                "vegan",
              
              # Data Wrangling 
                "tidyverse",
                "forcats",
                "readxl",
                "lubridate",
                "parallel",
                "broom",
                "recipes",
              
              # Data Viz
                "flextable",
                "ggplot2",
                "ggbeeswarm",
                "ggExtra",
                "ggrepel",
                "gridExtra",
                "RColorBrewer",
              
              # Statistics
                "MASS",
                "rcompanion",
                "glmmTMB",
                "lme4",
                "car",
                "nortest",
                "zCompositions",
                "caret"
                
              )) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}


## NonBioconductor Packages

install.packages("easystats", repos = "https://easystats.r-universe.dev")

### Additional packages to facilitate easystats, managed by BioManager
for (pkg in c('RUnit', 'wk', 'TTR', 'misc3d', 'mi', 'kutils', 'corrplot', 
              'lazyeval', 'elliptic', 'contfrac', 'deSolve', 'dygraphs', 
              'markdown', 'shinythemes', 'threejs', 'xts', 'svGUI', 
              'classInt', 's2', 'units', 'ECOSolveR', 'scs', 'osqp', 
              'operator.tools', 'furrr', 'qvcalc', 'relimp', 'flock', 
              'RhpcBLASctl', 'quantmod', 'FNN', 'multicool', 'plot3D', 
              'pracma', 'fdrtool', 'sem', 'XML', 'lisrelToR', 'rockchalk', 
              'statnet.common', 'ggsci', 'cowplot', 'ggsignif', 'polynom', 
              'rstatix', 'rpf', 'crosstalk', 'polyclip', 'hypergeo', 'flexmix',
              'mnormt', 'tmvnsim', 'StanHeaders', 'shinystan', 'nleqslv', 
              'maxLik', 'glmmML', 'miscTools', 'coneproj', 'svDialogs', 
              'scoringRules', 'estimability', 'BiasedUrn', 'pander', 
              'dreamerr', 'clue', 'CVXR', 'DEoptim', 'pbmcapply', 'Rcsdp', 
              'RSpectra', 'gamlss.data', 'gamlss.dist', 'plotrix', 'msm', 'sf',
              'pbivnorm', 'formula.tools', 'corpcor', 'cubature', 'arm',
              'broom.mixed', 'LaplacesDemon', 'dfidx', 'prediction', 
              'margins', 'gnm', 'ucminf', 'jtools', 'collapse', 'polspline',
              'fastGHQuad', 'mitools', 'lpSolve', 'bigassertr', 'bigparallelr',
              'nabor', 'fracdiff', 'tseries', 'urca', 'moments', 'metadat', 
              'mathjaxr', 'diptest', 'ks', 'enrichwith', 'RLRsim', 'oompaBase', 
              'oompaData', 'kernlab', 'CVST', 'qgraph', 'semPlot', 'glasso', 
              'network', 'GGally', 'fitdistrplus', 'ggpubr', 'matrixcalc', 
              'OpenMx', 'dendextend', 'DT', 'ellipse', 'flashClust', 'leaps', 
              'scatterplot3d', 'prabclus', 'SparseGrid', 'changepoint', 'cpm', 
              'gmp', 'Rmpfr', 'SuppDists', 'kSamples', 'BWStest',
              'distributional', 'fmsb', 'sjlabelled', 'sjmisc', 'reshape',
              'mc2d', 'HDInterval', 'ggforce', 'graphlayouts', 'Brobdingnag', 
              'inline', 'admisc', 'tweenr', 'evd'
              )
     ) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}


## Other Microbiome Packages -----------------------------------------------

devtools::install_github("ggloor/CoDaSeq/CoDaSeq")
devtools::install_github("kstagaman/phyloseqCompanion")