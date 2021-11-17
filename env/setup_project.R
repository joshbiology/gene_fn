#setup_project
#Goal: Setup necessary packages and validate folder structure to hold output figures.
#Follow each step to validate your environment.


# Step 1: Check installs --------------------------------------------------


### Check for and install missing packages

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

if (length(setdiff(pkgs, rownames(installed_pkgs)) > 0)) {
  install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

if (!("FastImputation" %in% rownames(installed_pkgs))) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/FastImputation/FastImputation_2.0.tar.gz", repos = NULL, type="source")
}

#This package needs to be installed from a specific historical link in CRAN


### Install bioconductor packages.
#Right now, I've removed any bioc dependencies in our repo for simplicity of deployment.
#If we need these packages eventually I'll turn this flag on.
check_bioc <- F

if (check_bioc) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(version = "3.14")
  
  
  bioc_packages <- c("AnnotationDbi", 
                     "Biobase", 
                     "BiocGenerics", 
                     "biomaRt", 
                     "fgsea", 
                     "GenomeInfoDb", 
                     "GenomicFeatures", 
                     "GenomicRanges", 
                     "IRanges", 
                     "mygene", 
                     "S4Vectors")
  
  if (length(setdiff(bioc_packages, rownames(installed_pkgs)) > 0)) {
    BiocManager::install(pkgs = setdiff(bioc_packages, installed_pkgs))
  }
}

# Step 2: Gather data --------------------------------------------------

#At this stage, you can either
#1. download the data manually:
#https://figshare.com/articles/dataset/Data_for_reproducing_Webster_figures/14960006
#and putting it into the gene_fn/data directory.

#2. Or, you can run gene_fn/load_data.sh, which automates #1 above.



# Step 3: Setup folders --------------------------------------------------
source("./env/setup_folders.R")


# Step 4: Test ProjectTemplate --------------------------------------------------

### Test load the package

library(ProjectTemplate)
load.project()

#The load.project command should run without errors. It may ask you to perform migrate.project() if it is missing any folders etc.


# Step 5: Setup graphics --------------------------------------------------
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




# Step 6: Run! ------------------------------------------------------------

#After these setup steps, you are ready to run items in the gene_fn/src folder.
#They are ordered to be run one after the other.
#supplement.R and portal_assets.R are helper scripts I wrote for 
#generating supplemental assets and file for our online portal. Users
#can skip those, as all those outputs are available here:
#https://figshare.com/articles/dataset/Webster_Supplemental_Output/14963561


