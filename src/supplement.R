#export for supplement
#Contains all export for supplemental tables in the Figshare repo.

library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "supplement")
create_output_folder(out_path)




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


# Function annot  --------------------------------------------------------------

#gprofiler
load("./cache/factor_gost_results.RData")

essential_genes <- read_tsv("./data/raw/depmap/public-19q4_v23-common-essentials.tsv")$gene %>% 
  convert_genes("cds_id", "symbol")

factor_genes <- map(1:webster_depmap$rank, function(x){
  loading <- get_gene_mat(webster_depmap)[,x]
  loading <- loading[order(desc(abs(loading)))]
  head(loading[loading != 0], 30) %>% names() %>% convert_cds_to_symbol()
})

frac_essential <- map_dbl(factor_genes, function(x) {
  length(intersect(x, essential_genes))/length(x)
})


factor_gost_results %>% 
  set_names( paste("V", 1:webster_depmap$rank, sep = "")) %>% 
  map_dfr(~.$result %>% arrange(p_value), .id = "Name") %>% 
  select(-query) %>% 
  left_join(tibble(Name = paste("V", 1:webster_depmap$rank, sep = ""),
         Frac_Essential = frac_essential,
         Loaded_Genes = factor_genes %>% map_chr(paste, collapse = " "))) %>% 
  select(Name, Frac_Essential, Loaded_Genes, everything()) %>% 
  select(-parents) %>% 
  write_tsv(file.path(out_path, "depmap_fn_annot_gprofiler.tsv"))

#Biomarkers
source("./munge/biomarker.R")
fn_models %>% 
  rename(Name = gene) %>% 
  write_tsv(file.path(out_path, "depmap_fn_biomarkers.tsv"))

# Other depmap resources --------------------------------------------------

#subcell enrichment scores
source("./munge/subcell.R")

load("./cache/max_nmf_df.RData")
write_tsv(max_nmf_df, file.path(out_path, "depmap_fn_subcell.tsv"))

load("./cache/nmf_projected.RData")

nmf_projected %>% 
  set_colnames(bioid_meta$Location) %>% 
  as_tibble(rownames = "Name") %>% 
  write_tsv(file.path(out_path, "depmap_fn_subcell_raw_matrix.tsv"))


#manual annotations
#See TableS3, Sheet 7 from Supplemental Information.


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




