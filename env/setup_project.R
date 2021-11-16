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

### Bioc
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("AnnotationDbi", 
                       "Biobase", 
                       "BiocGenerics", 
                       "biomaRt", 
                       "fgsea", 
                       "GenomeInfoDb", 
                       "GenomicFeatures", 
                       "GenomicRanges", 
                       "IRanges", 
                       "mygene", 
                       "S4Vectors"))

### Gather data
#At this stage, you can either
#1. download the data manually:
#https://figshare.com/articles/dataset/Data_for_reproducing_Webster_figures/14960006
#and putting it into a directory called gene_fn/data

#2. You can run gene_fn/load_data.sh, which is used by the Dockerfile to do the same as #1 above.


### Test load the package

library(ProjectTemplate)
load.project()

#The load.project command should run without errors. It may ask you to perform migrate.project() if it is missing any folders etc.



### Setup folders
source("./env/setup_folders.R")
