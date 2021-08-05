library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "05-function_annot")



# Output gene fn heatmaps -------------------------------------------------

source("./munge/webster_depmap.R")
source("./munge/depmap_sample_info.R")


plot_gene_loadings <- function(index, max_loading = F, flip = F) {
  gene_df <- get_gene_mat(webster_depmap)[,index] %>% 
    enframe("CDS_ID", "Loadings") %>% 
    mutate(symbol = convert_cds_to_symbol(CDS_ID))
  
  
  if (flip) {
    direction_indicator  <-  -1
    
  }
  
  else
    direction_indicator  <-  1
  
  gene_df <- gene_df %>% 
    mutate(Loadings = Loadings * direction_indicator,
           Rank = order(order(Loadings, decreasing = T))) %>% 
    arrange(desc(Loadings)) %>% 
    filter(Rank <= 10)
  
  loading_mat <- get_gene_mat(webster_depmap)[gene_df$CDS_ID,index] %>% matrix(ncol = 1) %>% 
    set_rownames(gene_df$symbol)
  
  loading_mat <- loading_mat * direction_indicator 
  
  if (max_loading == F) {
    max_loading <- max(abs(loading_mat))
  }
  
  g2 <- pheatmap::pheatmap(loading_mat, cluster_cols = F, cluster_rows = F, show_colnames = F, 
                           filename =  file.path(out_path, paste("gene_loading_", index, ".pdf", sep = "")),
                           cellwidth =20, cellheight = 20, color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), 
                           breaks = seq(-max_loading, max_loading, length.out=101))
}

walk(1:webster_depmap$rank, plot_gene_loadings)

plot_cell_loadings <- function(index, max_loading = F, annot = NA) {
  
  cell_df <- get_cell_mat(webster_depmap)[,index] %>% 
    enframe("DepMap_ID", "Loadings") %>% 
    arrange(Loadings) %>% 
    mutate(Name = convert_depmap_id_to_ccle(DepMap_ID),
           Rank = order(order(Loadings, decreasing = F))) %>% 
    filter(Rank <= 10) %>% 
    arrange((Loadings))
  
  loading_mat <- get_cell_mat(webster_depmap)[cell_df$DepMap_ID,index] %>% matrix(ncol = 1) %>% 
    set_rownames(cell_df$Name)
  
  
  if (!max_loading) {
    max_loading <- max(abs(loading_mat))
  }
  
  if(is.na(annot)) {
    
    
    
    g1 <- pheatmap::pheatmap(loading_mat, cluster_cols = F, cluster_rows = F, show_colnames = F, 
                             filename =  file.path(out_path, paste("cell_loading_", index, ".pdf", sep = "")),
                             cellwidth =20, cellheight = 20, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(),
                             breaks = seq(-max_loading, max_loading, length.out=101))
  }
  
  else
    g1 <- pheatmap::pheatmap(loading_mat, cluster_cols = F, cluster_rows = F, show_colnames = F, annotation_row = annot,
                             filename =  file.path(out_path, paste("cell_loading_", index, ".pdf", sep = "")),
                             cellwidth =20, cellheight = 20, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(),
                             breaks = seq(-max_loading, max_loading, length.out=101))
  
  return(loading_mat)
}

walk(1:webster_depmap$rank, plot_cell_loadings)


# Gene fn annotation  -------------------------------------------------
essential_genes <- read_tsv("./data/raw/depmap/public-19q4_v23-common-essentials.tsv")$gene %>% 
  convert_genes("cds_id", "symbol")

hart_essentials <- read_tsv("./data/raw/hart/")

gene_loadings_df <- webster_depmap %>% 
  get_gene_mat() %>% 
  set_colnames(paste("V", 1:webster_depmap$rank, sep = "")) %>% 
  as_tibble(rownames = "Name") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V"))


library(gprofiler2)



factor_genes <- map(1:webster_depmap$rank, function(x){
  loading <- get_gene_mat(webster_depmap)[,x]
  loading <- loading[order(desc(abs(loading)))]
  head(loading[loading != 0], 30) %>% names() %>% convert_cds_to_symbol()
})

frac_essential <- map_dbl(factor_genes, function(x) {
  length(intersect(x, essential_genes))/length(x)
})

factor_gost_results <- map(factor_genes, function(x) {
  gprofiler2::gost(query = c(query = x, organism = "hsapiens", ordered_query = T))
})

ProjectTemplate::cache("factor_gost_results")

# Aggregate and name ------------------------------------------------------\
require(gprofiler2)
load("./cache/factor_gost_results.RData")

#Plot gost output
#https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html

annotate_fn_plot <- function(index, num_terms = 10) {
  gost <- factor_gost_results[[index]]
  
  file_out <- file.path(out_path, paste("gprofiler_", index, ".pdf", sep = ""))
  top_terms <- gost$result %>% arrange(p_value) %>% slice(1:num_terms)
  publish_gostplot(gostplot(gost, interactive = FALSE), highlight_terms = top_terms, width = 10, height = 8, filename = file_out)
}

walk(1:webster_depmap$rank, annotate_fn_plot)


#Export top enrichment
factor_gost_results %>% 
  map_dfr(~.$result %>% arrange(p_value) %>% slice(1)) %>% 
  mutate(Name = paste("V", 1:webster_depmap$rank, sep = ""),
         Frac_Essential = frac_essential,
         Loaded_Genes = factor_genes %>% map_chr(paste, collapse = " "))%>% 
  select(Name, Frac_Essential, Loaded_Genes,everything()) %>% 
  select(-parents) %>% 
  write_tsv(file.path(out_path, "depmap_fn_annot_gprofiler.tsv"))



# STRINGdb ----------------------------------------------------------------
#https://bioconductor.org/packages/devel/bioc/manuals/STRINGdb/man/STRINGdb.pdf




# Number of pleiotropic genes ---------------------------------------------

threshold <- 0.25
loadings_df <- webster_depmap %>%
  get_gene_mat_stdev() %>% 
  as_tibble(rownames=  "Gene") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", -Gene) %>% 
  filter(abs(Loading) > threshold) %>% 
  arrange(Gene, Function) 

pleiop_gene_count <-loadings_df %>% count(Gene) %>% 
  left_join(gene_recon_webster)


(pleiop_gene_count %>% filter(n != 1) %>% nrow())/(nrow(pleiop_gene_count))

factor_pleiotrop_df <- loadings_df %>% 
  left_join(loadings_df, by = "Gene") %>% 
  filter(Function.x != Function.y, Function.x < Function.y) %>% 
  left_join(gene_recon_webster) %>% 
  filter(Recon_Pearson > 0.4)

nrow(factor_pleiotrop_df)

plot_pleiotropy_network <- function(fns) {
  require(igraph)
  require(tidygraph)
  require(ggraph)
  tmp_graph <- as_tbl_graph(igraph::graph_from_data_frame(factor_pleiotrop_df %>% 
                                                            filter(Function.x %in% fns & Function.y%in% fns) %>%
                                                            count(Function.x, Function.y) %>% 
                                                            arrange(desc(n)), directed = F))
  
  ggraph(tmp_graph, layout = "circle") + 
    geom_edge_link(color = "gray90", alpha = 1, aes(edge_width = n)) + 
    geom_node_text(aes(label = name), size = 5, color = "black") +
    geom_node_point(size = 2.5, color = "black")
  
}

#Shoc2 function selection
factor_pleiotrop_df %>% select(Function.x, Function.y) %>% 
  count(Function.x, Function.y) %>% 
  filter(Function.x %in% c("V138", "V209", "V162", "V8") | Function.y %in% c("V138", "V209", "V162", "V8")) %>% View()

#SHOC2
plot_pleiotropy_network(c("V138", "V209", "V162", "V8", "V215", "V6", "V15")) +
  ggsave(file.path(out_path, "shoc2_pleiop.pdf"), width = 4, height = 3.5, device = cairo_pdf)

#Sterol
plot_pleiotropy_network(c("V51", "V30")) +
  ggsave(file.path(out_path, "sterol_pleiop.pdf"), width = 4, height = 3.5, device = cairo_pdf)

#retro
plot_pleiotropy_network(c("V17", "V169", "V11")) +
  ggsave(file.path(out_path, "retro_pleiop.pdf"), width = 4, height = 3.5, device = cairo_pdf)

#scar
plot_pleiotropy_network(c("V52", "V6", "V34")) +
  ggsave(file.path(out_path, "scar_pleiop.pdf"), width = 4, height = 3.5, device = cairo_pdf)


# Biomarker ---------------------------------------------------------------
#Plot bargraphs of Omics associations.


source("./munge/biomarker.R")

plot_biomarkers <- function(index, trim_length = 55) {
  
  biom_df <- fn_models %>% 
    filter(gene == paste("V", index, sep = "")) %>% 
    mutate(Feature_Name = str_sub(Feature_Name, 1, trim_length)) %>% 
    mutate(Feature_Name = tidytext::reorder_within(Feature_Name, importance, model))
  
  labels <- biom_df %>% 
    group_by(model) %>% 
    summarize(Label = paste("Pearson =", max(pearson) %>% round(digits = 2)))
    
  biom_df %>% 
    ggplot(aes(Feature_Name, importance)) +
    geom_col(aes(fill = Type)) +
    tidytext::scale_x_reordered() +
    coord_flip() +
    facet_wrap(~model, scales = "free") +
    geom_text(data = labels, aes(label = Label, x = -Inf, y = -Inf), vjust = -1, hjust=-1) +
    ggtitle(paste("V", index, sep = "")) +
    ggsave(file.path(out_path, paste("biomarker_", index, "_", "scores", ".pdf", sep = "")), width = 15, height = 8, device= cairo_pdf)
}

walk(1:webster_depmap$rank, plot_biomarkers)
