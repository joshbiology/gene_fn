# Gene Function (Webster)
This is the repository for reproducing figures and analyses from:
"Sparse dictionary learning recovers pleiotropy from human cell fitness screens", by Pan et al. 

This R codebase uses pre-factorized MATLAB objects as the basis for many figures. If you would like try the factorization yourself, that MATLAB codebase is here: https://github.com/joshbiology/graph_dictionary_learning.

# Setup

We include a table of the required packages in order to run all of the scripts included. If these packages are missing from your R environment, they will flag missing package errors during runtime. One can either install all the packages up front, or run individual scripts of interest and manage package dependencies on the fly.

We are also preparing a docker image of the frozen R and package dependnecies that one can use as the basis for running these scripts. 

# Obtaining data

The figshare repository [https://doi.org/10.6084/m9.figshare.14960006] consists of a single zipped folder, named 'data', that one should download and copy (or mount) into the root folder. The structure of the folder will look as follows:

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

# Figure to directory mapping
```
Figure 1 | 01-synthetic_example
Figure 2 | 02-genotoxic
Figure 3 | 03-depmap_preprocess, 04-depmap_postprocess, 05-function_annot, 06-sparse_codes, 08-landscape
Figure 4 | 07-complexes
Figure 5 | 08-landscape
Figure 6 | 09-projection
```

