library(ProjectTemplate); load.project()
source("./munge/webster_depmap.R")

out_path <- file.path(".", "output", "07-complexes")
create_output_folder(out_path)

# Functions ----------------------------------------------------------------




# Import data -------------------------------------------------------------


genesets <- list.files(file.path(".", "data", "raw", "genesets"), full.names = T) %>% 
  map(~read_tsv(.) %>% mutate(entrezgene = convert_genes_mygeneinfo(Gene ,"entrezgene"),
                              CDS_ID = convert_genes(entrezgene, "entrez_id", "cds_id"))) %>% 
  set_names(list.files(file.path(".", "data", "raw", "genesets"), full.names = F) %>% str_sub(1, -5)) 






# Plot complex heatmaps ---------------------------------------------------

#Identify factors
plot_loading_heatmap <- function(df, seed_gene, num_of_functions, out_name, loading_cap =0.75) {
  focus_df <- df %>% select(CDS_ID, Geneset) %>% arrange(Geneset)
  
  selected_functions <- get_gene_loadings(webster_depmap, convert_genes(seed_gene, "symbol", "cds_id")) %>% 
    arrange(-Loading) %>% 
    slice(1:num_of_functions) %>% 
    pull(Function)
  
  genes_tmp <- intersect(focus_df$CDS_ID, avana_19q4_webster %>% colnames)
  
  loading_mat <- webster_depmap %>% get_gene_mat_stdev() %>% magrittr::extract(genes_tmp,selected_functions) 
  
  pheatmap::pheatmap(loading_mat,cluster_cols = F, cluster_rows = F, 
                     annotation_row = focus_df %>% column_to_rownames("CDS_ID"), breaks = seq(0, loading_cap, length.out=101),
                     color = colorRampPalette(c("white", "#5ecbf5"))(101), cellwidth = 15, cellheight = 15, filename = out_name)
}




plot_loading_heatmap(genesets$staga_atac, "KAT2A", 2, file.path(out_path, "staga.pdf"))
plot_loading_heatmap(genesets$baf, "SMARCA4", 3, file.path(out_path, "baf.pdf"))
plot_loading_heatmap(genesets$mediator, "MED1", 2, file.path(out_path, "med.pdf"))
plot_loading_heatmap(genesets$ints, "C7orf26", 3, file.path(out_path, "ints.pdf"))

