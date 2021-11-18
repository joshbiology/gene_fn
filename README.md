# Gene Function (Webster)
This is the repository for reproducing figures and analyses from:
"Sparse dictionary learning recovers pleiotropy from human cell fitness screens", by Pan et al. 2021, available on [arXiv](https://arxiv.org/abs/2111.06247).

This R codebase uses pre-factorized MATLAB objects as the basis for many figures. If you would like try the factorization yourself, that MATLAB codebase is here: https://github.com/joshbiology/graph_dictionary_learning.

# Setup

## Clone this repository

From the command line:

```
$ git clone https://github.com/joshbiology/gene_fn.git
```


## Acquire data

You can download the data manually at
https://figshare.com/articles/dataset/Data_for_reproducing_Webster_figures/14960006
and place it into a gene_fn/data folder.

Or, you can run gene_fn/load_data.sh, which automates the above.


The figshare repository consists of a single zipped folder, named 'data' structured as follows:

```
data
├── interim
│   ├── biomarker
│   └── matlab
│   	├── depmap_deep
│   	├── depmap_denoise
│   	├── depmap_grid
│   	├── depmap_wide
│   	├── durocher
│   	├── durocher_no_graph
│   	└── synthetic
└── raw
    ├── depmap
    ├── durocher
    ├── genesets
    ├── hart
    ├── humancellmap
    ├── nusiance_genesets
    └── prism
```

If you are interested in only the (1) input data matrices for Webster's factorizations and/or (2) The Webster output that were the basis of the figures, we created a second repository that holds these supplemental tables [https://doi.org/10.6084/m9.figshare.14963561].


## Acquire software dependencies

The code has been tested on R version 4.1.2. 

There are two ways to acquire needed packages: install locally, or use our Docker image.

If you are running locally, skip to Setup Project below.

If you wish to use the docker image, install [Docker](https://www.docker.com) and run:

```
$ docker pull jlpan2/gene_fn
```

Then, you can mount the cloned directory and start a docker container:

```
docker run -v $[PATH TO CLONED DIR]:/gene_fn -it jlpan2/gene_fn

```

## Setup project

This script contains everything you need to set up your project locally. In an R session
in local clone of gene_fn, run 

```
source("./env/setup_project.R")
```
And follow through the instructions in the comments. You may need to perform
local operations to enable graphics. If you are using a docker container, you should
still run this step but those dependencies should be available to you.

## Run the src directory.
You can run each file in the src directory, or run them all using:

```
ProjectTemplate::run.project()
```

If using Docker, make sure to allocate plenty of memory for your containers, as 
processes will crash once memory limits are reached.

# Figure to directory mapping
```
Figure 1 | 01-synthetic_example
Figure 2 | 02-genotoxic
Figure 3 | 03-depmap_preprocess, 04-depmap_postprocess, 05-function_annot, 06-sparse_codes, 08-landscape
Figure 4 | 07-complexes
Figure 5 | 08-landscape
Figure 6 | 09-projection
```


# sessionInfo()

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_0.5.1          stringr_1.4.0          dplyr_1.0.7            purrr_0.3.4            readr_2.1.0            tidyr_1.1.4            ggplot2_3.3.5          tidyverse_1.3.1       
 [9] umap_0.2.7.0           ff_4.0.5               bit_4.0.4              magrittr_2.0.1         readxl_1.3.1           plyr_1.8.6             reshape2_1.4.4         ProjectTemplate_0.10.2
[17] tibble_3.1.6           digest_0.6.28         

loaded via a namespace (and not attached):
 [1] matrixStats_0.61.0  fs_1.5.0            lubridate_1.8.0     bit64_4.0.5         RColorBrewer_1.1-2  httr_1.4.2          tools_4.1.2         backports_1.3.0     utf8_1.2.2         
[10] R6_2.5.1            DBI_1.1.1           lazyeval_0.2.2      colorspace_2.0-2    withr_2.4.2         tidyselect_1.1.1    gridExtra_2.3       tictoc_1.0.1        compiler_4.1.2     
[19] cli_3.1.0           rvest_1.0.2         xml2_1.3.2          plotly_4.10.0       scales_1.1.1        askpass_1.1         R.utils_2.11.0      pkgconfig_2.0.3     htmltools_0.5.2    
[28] parallelly_1.28.1   dbplyr_2.1.1        fastmap_1.1.0       htmlwidgets_1.5.4   rlang_0.4.12        rstudioapi_0.13     generics_0.1.1      jsonlite_1.7.2      vroom_1.5.6        
[37] R.oo_1.24.0         R.matlab_3.6.2      Matrix_1.3-4        Rcpp_1.0.7          munsell_0.5.0       fansi_0.5.0         reticulate_1.22     viridis_0.6.2       lifecycle_1.0.1    
[46] R.methodsS3_1.8.1   furrr_0.2.3         stringi_1.7.5       grid_4.1.2          parallel_4.1.2      listenv_0.8.0       crayon_1.4.2        lattice_0.20-45     haven_2.4.3        
[55] hms_1.1.1           pillar_1.6.4        igraph_1.2.8        codetools_0.2-18    reprex_2.0.1        glue_1.5.0          data.table_1.14.2   BiocManager_1.30.16 modelr_0.1.8       
[64] png_0.1-7           vctrs_0.3.8         tzdb_0.2.0          cellranger_1.1.0    gtable_0.3.0        openssl_1.4.5       future_1.23.0       assertthat_0.2.1    broom_0.7.10       
[73] RSpectra_0.16-0     viridisLite_0.4.0   pheatmap_1.0.12     mltools_0.3.5       globals_0.14.0      ellipsis_0.3.2     

```

