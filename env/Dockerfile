FROM rocker/r-ver:4.1.1
LABEL maintainer="joshpan"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y git-core \
	libcurl4-openssl-dev \
	libglpk-dev \
	libgmp-dev \
	libssl-dev \
	libxml2-dev \
	make \
	pandoc \
	pandoc-citeproc \
  libgsl-dev

RUN ["install2.r", "bitops", "caret", "coop", "coRanking", "doParallel", "dplyr", "entropy", "fastICA", "FastImputation", "forcats", "foreach", "furrr", "future", "ggExtra", "ggplot2", "ggraph", "gprofiler2", "httr", "igraph", "iterators", "lattice", "magrittr", "Matrix", "matrixStats", "Metrics", "mltools", "modes", "pheatmap", "plyr", "ProjectTemplate", "purrr", "R.matlab", "RColorBrewer", "Rcpp", "RcppZiggurat", "RCurl", "readr", "readxl", "remotes", "Rfast", "rfigshare", "sn", "stringr", "tibble", "tictoc", "tidygraph", "tidyr", "tidytext", "tidyverse", "umap", "viridis", "viridisLite", "yardstick"]

RUN ["install2.r", "-r https://bioconductor.org/packages/3.10/bioc -r https://bioconductor.org/packages/3.10/data/annotation -r https://bioconductor.org/packages/3.10/data/experiment -r https://bioconductor.org/packages/3.10/workflows", "AnnotationDbi", "Biobase", "BiocGenerics", "biomaRt", "fgsea", "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "IRanges", "mygene", "S4Vectors"]

RUN ["installGithub.r", "r-lib/scales@a33f8ef9d2930f4940f9c694ba6ba79480965cec", "ropensci/plotly@a689263996ae0938c2f18987877924178dba9cf6", "rstudio/renv@48b65fc927fa975971e480dc58f44d86f2243208"]

RUN apt-get install libxtst6 -y

WORKDIR /gene_fn/
