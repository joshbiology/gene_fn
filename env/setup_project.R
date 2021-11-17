#setup_project
#Goal: For users not using Docker, setup necessary packages locally.

###Manual install
#This package needs to be installed from a specific historical link in CRAN
install.packages("https://cran.r-project.org/src/contrib/Archive/FastImputation/FastImputation_2.0.tar.gz", repos = NULL, type="source")

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

#2. You can run gene_fn/load_data.sh, which automates #1 above.


### Setup folders
source("./env/setup_folders.R")

### Test load the package

library(ProjectTemplate)
load.project()

#The load.project command should run without errors. It may ask you to perform migrate.project() if it is missing any folders etc.



### Setup graphics
#https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
#Many plots use cairo_pdf to export rich PDF graphics.
#There are different setup requirements on Windows, R and Linux.
#If you are seeing plots not show up, and errors that look like this, you might need to install X11 (Linux or Mac)

#Example warning 1:
#In grSoftVersion() :
#  unable to load shared object '/usr/local/lib/R/modules//R_X11.so':
#  libXt.so.6: cannot open shared object file: No such file or directory

#Example warning 2:
#In (function (filename = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf",  ... :
#               failed to load cairo DLL
              
#We have taken care of this in our Docker dependencies, so if you're using Docker this should be fine.
