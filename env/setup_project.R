#setup_project
#Goal: For users not using Docker, setup necessary packages locally.
### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("bitops", 
           "caret", 
           "coop", 
           "coRanking", 
           "doParallel", 
           "dplyr", 
           "entropy", 
           "fastICA", 
           "FastImputation", 
           "forcats", 
           "foreach", 
           "furrr", 
           "future", 
           "ggExtra", 
           "ggplot2", 
           "ggraph", 
           "gprofiler2", 
           "httr", 
           "igraph", 
           "iterators", 
           "lattice", 
           "magrittr", 
           "Matrix", 
           "matrixStats", 
           "Metrics", 
           "mltools", 
           "modes", 
           "pheatmap", 
           "plyr", 
           "ProjectTemplate", 
           "purrr", 
           "R.matlab", 
           "RColorBrewer", 
           "Rcpp", 
           "RcppZiggurat", 
           "RCurl", 
           "readr", 
           "readxl", 
           "remotes", 
           "Rfast", 
           "rfigshare", 
           "sn", 
           "stringr", 
           "tibble", 
           "tictoc", 
           "tidygraph", 
           "tidyr", 
           "tidytext", 
           "tidyverse", 
           "umap", 
           "viridis", 
           "viridisLite", 
           "yardstick"
)

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
  install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

### Test load the package

library(ProjectTemplate)
load.project()
