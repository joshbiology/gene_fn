#for portal export. For each dataset XXX, assets need to be:
#1. XXX_umap.tsv - UMAP coordinates. Name, X, Y, Type (Gene or Function).
#2. XXX_gene_to_function.tsv - Loadings. Name, Column per Function, with gene loadings inside.
#3. XXX_gene_meta.tsv - Info for tooltips and gene cards. Loading_1, Function_1, etc. for gene cards, and all the relevant info for tooltips.
#4. XXX_fn_meta.tsv - Info for tooltips and fn cards.
#5. XXX_collections.tsv - Hierarchy for collections in fn selector. Function, Subcollection and Collection.
#6. XXX_svg_FUNCTIONNAME.svg - Small multiple of UMAP for fn selector.

#For DepMap specifically:
#7. Systematic ID to Manual name - key value pair, Name, Display_Name
#8. Depmap matrices - for biomarker pipeline.

# Depmap ------------------------------------------------------------------


library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "portal_assets")


#Depmap
load("./cache/avana_19q4_webster.RData")
source("./munge/webster_depmap.R")


# Depmap UMAP -------------------------------------------------------------
load('./cache/plotting_df.RData')

plotting_df %>% 
  select(-Gene_Name) %>% 
  write_tsv(file.path(out_path, "depmap_umap.tsv"))

# Depmap gene to function -------------------------------------------------

gene_loadings <- webster_depmap %>% 
  get_gene_mat() %>% 
  set_colnames(paste("V", 1:webster_depmap$rank, sep = "")) %>% 
  as_tibble(rownames = "Name") 

write_tsv(gene_loadings, file.path(out_path, "depmap_gene_to_function.tsv"))


# Depmap small multiples --------------------------------------------------
gene_loadings_df <- gene_loadings %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V"))

walk(1:webster_depmap$rank, function(x){
  
  factor_name = paste("V", x, sep = "")
  
  tmp_df <-   plotting_df %>% 
    left_join(gene_loadings_df) %>% 
    filter(Factor == factor_name | Name == factor_name) 
  
  tmp_df %>% 
    ggplot(aes(X, Y, color = Loading, shape = Type)) +
    geom_point(data = subset(tmp_df, Type == "Gene" & Loading == 0),alpha = 0.85) +
    geom_point(data = subset(tmp_df, Type == "Gene" & Loading != 0),alpha = 0.85) +
    geom_point(data = subset(tmp_df, Type == "Function"),alpha = 0.85) +
    scale_shape_manual(values=c(17, 16)) +
    scale_color_gradient2(low = "#BD6C33", mid = "gray90", high =  "#00AEEF", midpoint = 0, limits=c(-10, 10), oob = scales::squish) +
    theme_void() +
    theme(legend.position = "NA") +
    ggsave(filename = sprintf("./output/portal_assets/depmap_svg_%s.svg", factor_name), width = 5, height = 5)
  
})

# Depmap gene metadata -------------------------------------------
# Subcellular localization info
source("./munge/subcell.R")

bioid_card <- bioid_loc %>%
  mutate(Name = convert_genes(entrezgene,"entrez_id", "cds_id")) %>% 
  dplyr::select(Name, Location)


# Portal page URL
#http://depmap.org
#https://www.genecards.org/Guide/AboutGeneCards
#https://www.ncbi.nlm.nih.gov/gene/
depmap_gene_urls <- tibble(Name = avana_19q4_webster %>% colnames(),
                    DepMap_URL = map_chr(Name, function(x) sprintf("https://depmap.org/portal/gene/%s?tab=overview", convert_cds_to_symbol(x))),
                    GeneCard_URL = map_chr(Name, function(x) sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", convert_cds_to_symbol(x))),
                    NIH_Gene_URL = map_chr(Name, function(x) sprintf("https://www.ncbi.nlm.nih.gov/gene/?term=%s", convert_cds_to_entrez(x)))) 


# Reconstruction score
gene_recon_webster

#Loadings
depmap_loading_fn_card <- gene_loadings_df %>% 
  filter(Loading != 0) %>% 
  group_by(Name) %>% 
  arrange(Name, desc(abs(Loading))) %>% 
  ungroup() %>% 
  mutate(Rank = rep(1:4, ncol(avana_19q4_webster))) %>% 
  pivot_wider(names_from = "Rank", values_from = c("Function", "Loading"))

#Combine
depmap_gene_meta <- depmap_loading_fn_card %>% 
  left_join(gene_recon_webster %>% rename("Name" = "Gene")) %>% 
  left_join(bioid_card) %>% 
  left_join(depmap_gene_urls)

write_tsv(depmap_gene_meta, file.path(out_path, "depmap_gene_meta.tsv"))

#consider adding PPI to the cards / tooltips
#https://bioconductor.org/packages/release/data/experiment/vignettes/simpIntLists/inst/doc/simpIntLists.R



# Depmap display name -----------------------------------------------------
load("./cache/factor_gost_results.RData")

depmap_gost_meta <- factor_gost_results %>% 
  map_dfr(function(x) x$result %>% arrange(p_value) %>% slice(1)) %>% 
  mutate(Name = paste("V", 1:webster_depmap$rank, sep = ""))

depmap_fn_display_name <- depmap_gost_meta %>% 
  select(Name, Display_Name = term_name) %>% 
  mutate(Display_Name = str_sub(Display_Name,1, 30)) %>% 
  as_tibble()

write_tsv(depmap_fn_display_name, file.path(out_path, "depmap_fn_display_name.tsv"))

# Depmap collections ----------------------------------------------------

load("./cache/max_nmf_df.RData")
load("./cache/bioid_colors.RData")


depmap_loc_meta <- max_nmf_df %>% 
  left_join(bioid_colors) %>% 
  group_by(Compartment) %>% 
  arrange(Compartment,desc(Specificity), Location) %>% 
  filter(Specificity > 0.33) %>% 
  slice(1:15) %>% 
  ungroup() %>% 
  


depmap_loc_collection <-depmap_loc_meta %>%
  arrange(Compartment, Specificity) %>% 
  mutate(Collection = "Subcellular") %>% 
  select(Name, Collection, Subcollection = Compartment)


#depmap_cancer_collection 

write_tsv(depmap_loc_collection, file.path(out_path, "depmap_collections.tsv"))


# Depmap fn metadata --------------------------------------------------
#ingest biomarkers
fn_biomarkers <- read_tsv("./data/interim/biomarker/latent-representation_v4-latent-var-ensemble.tsv") %>% 
  select(1:8)  %>% filter(best) %>% 
  rename(Name = gene, Biomarker_Pearson = pearson) %>% 
  select(-model,-best)

depmap_fn_meta <- depmap_gost_meta %>% 
  select(Name, Term_Name = term_name, Term_ID = term_id, Term_P_Value = p_value) %>% 
  left_join(max_nmf_df %>%   left_join(bioid_colors) %>% 
    select(Name, Compartment, Location, Loc_Specificity = Specificity)) %>% 
  mutate(Compartment= factor(Compartment, levels = c("Mito", "Nucleus", "Membrane", "Trafficking", "ER" , "Cytoskeleton", "Misc.")))  %>% 
  left_join(fn_biomarkers)


write_tsv(depmap_fn_meta, file.path(out_path, "depmap_fn_meta.tsv"))

# DepMap matrices ---------------------------------------------------------
#Exported for biomarker analysis
avana_19q4_webster %>% 
  as_tibble(rownames = "Cell_Line") %>% 
  write_csv(file.path(out_path, paste(Sys.Date(), "avana_19q4_webster.csv", sep = "_")))


get_cell_mat(webster_depmap) %>% 
  as_tibble(rownames = "Cell_Line") %>% 
  write_csv(file.path(out_path, paste(Sys.Date(), "dictionary_matrix_19q4_webster.csv", sep = "_")))






# Durocher ----------------------------------------------------------------

# Durocher gene metadata ----------------------------------------------------------------

durocher_g2f <- read_tsv("./output/portal_assets/durocher_gene_to_function.tsv") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", -Name) %>% 
  filter(Loading != 0) %>% 
  group_by(Name) %>% 
  arrange(Name, desc(abs(Loading))) %>% 
  ungroup() %>% 
  mutate(Rank = rep(1:2, 304)) %>% 
  pivot_wider(names_from = "Rank", values_from = c("Function", "Loading"))

durocher_recon <- read_tsv("./output/02-genotoxic/durocher_recon_pearson.tsv")

olivieri_genes <- read_excel("./data/raw/durocher/cell/mmc4-2.xlsx", sheet = 1, skip = 1) %>%
  pivot_longer(names_to = "Literature_Pathway", values_to = "Name", everything()) %>% 
  group_by(Name) %>% 
  summarize(Literature_Pathway = paste(Literature_Pathway, collapse = "; "))

durocher_gene_urls <- tibble(Name = durocher_g2f %>% pull(Name),
                           GeneCard_URL = map_chr(Name, function(x) sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", x)))

durocher_gene_meta <- durocher_g2f %>% 
  left_join(durocher_recon) %>% 
  left_join(olivieri_genes) %>% 
  left_join(durocher_gene_urls)

write_tsv(durocher_gene_meta, file.path(out_path, "durocher_gene_meta.tsv"))

# Durocher fn metadata ----------------------------------------------------------------
#durocher_cell_mat <- $
  
durocher_fn_top <- read_tsv("./output/portal_assets/durocher_cell_mat.tsv") %>% 
  pivot_longer(names_to = "Treatment", values_to = "Dependency", -Name) %>% 
  group_by(Name) %>% 
  arrange(Name, desc(abs(Dependency))) %>% 
  slice(1:2) %>% 
  mutate(Rank = 1:2) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Rank", values_from = c("Treatment", "Dependency"))

write_tsv(durocher_fn_top, file.path(out_path, "durocher_fn_meta.tsv"))
