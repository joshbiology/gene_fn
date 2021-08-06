#export for supplement
#Contains all export for supplemental tables in the Figshare repo.

library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "supplement")




# Genotoxic ---------------------------------------------------------------

source("./munge/webster_genotoxic.R")


#Input data matrix
write_tsv(genotoxic_input %>% as_tibble(rownames = "Treatment"), file.path(out_path, "genotoxic_input.tsv"))


#Factorized output
#Gene loadings
get_gene_mat_stdev(webster_genotoxic) %>% 
  magrittr::extract(,function_order) %>% 
  set_colnames(function_names[function_order]) %>% 
  as_tibble(rownames = "Gene") %>% 
  write_tsv(file.path(out_path, "genotoxic_gene_loadings.tsv"))

#Dictionary matrix
get_cell_mat(webster_genotoxic) %>% 
  magrittr::extract(,function_order) %>% 
  set_colnames(function_names[function_order]) %>% 
  as_tibble(rownames = "Treatment") %>% 
  write_tsv(file.path(out_path, "genotoxic_dictionary.tsv"))

#The following are generated from portal_assets.R:
#Gene Annotations
#Function Annotations
#Embedding coordinates

# Cancer Dependency Map (DepMap) ------------------------------------------
source("./munge/webster_depmap.R")


#Input data matrix
write_tsv(avana_19q4_webster %>% as_tibble(rownames = "DepMap_ID"), file.path(out_path, "depmap_input.tsv"))



#Factorized output
#Gene loadings
get_gene_mat_stdev(webster_depmap) %>% 
  as_tibble(rownames = "Gene") %>% 
  write_tsv(file.path(out_path, "depmap_gene_loadings.tsv"))

#Dictionary matrix
get_cell_mat(webster_depmap) %>% 
  as_tibble(rownames = "DepMap_ID") %>% 
  write_tsv(file.path(out_path, "depmap_dictionary.tsv"))




#Cell Line annotations
source("./munge/depmap_sample_info.R")

sample_info %>% 
  filter(DepMap_ID %in% rownames(get_cell_mat(webster_depmap))) %>% 
  write_tsv(file.path(out_path, "depmap_cell_line_info.tsv"))




#The following are generated from portal_assets.R:
#Function Annotations
#Gene Annotations
#UMAP embedding coordinates


# Other depmap resources --------------------------------------------------

#subcell enrichment scores
load("./cache/max_nmf_df.RData")
write_tsv(max_nmf_df, file.path(out_path, "depmap_fn_subcell.tsv"))

#gprofiler enrichment
#Exported as part of 05-function_annot.R

#manual annotations

# PRISM annotations -------------------------------------------------------

#Combined primary embedding
load("./cache/impute_df.RData")
impute_df %>% 
  write_tsv(file.path(out_path, "prism_embedding.tsv"))


#Primary input data matrix
load("./cache/prism_primary_imputed.RData")
write_tsv(prism_primary_imputed %>% as_tibble(rownames = "DepMap_ID"), file.path(out_path, "prism_primary_imputed.tsv"))


#Compound loadings
load("./cache/prism_primary_omp.RData")
write_tsv(prism_primary_omp %>% as_tibble(rownames = "column_name"), file.path(out_path, "prism_primary_omp.tsv"))


#Compound meta
load("./cache/primary_compound_filtered_meta.RData")
write_tsv(primary_compound_filtered_meta, file.path(out_path, "prism_primary_meta.tsv"))

#Compound results
load("./cache/proj_results.RData")
write_tsv(proj_results, file.path(out_path, "prism_primary_proj_results.tsv"))

#Secondary input data matrix
load("./cache/prism_secondary_imputed.RData")
write_tsv(prism_secondary_imputed %>% as_tibble(rownames = "DepMap_ID"), file.path(out_path, "prism_secondary_imputed.tsv"))

#Secondary Compound loadings
load("./cache/prism_secondary_omp.RData")
write_tsv(prism_secondary_omp %>% as_tibble(rownames = "column_name"), file.path(out_path, "prism_secondary_omp.tsv"))

#Secondary compound meta
load("./cache/secondary_meta.RData")
write_tsv(secondary_meta, file.path(out_path, "prism_secondary_meta.tsv"))


#Compound results
load("./cache/proj_results_secondary.RData")
write_tsv(proj_results_secondary, file.path(out_path, "prism_secondary_proj_results.tsv"))




