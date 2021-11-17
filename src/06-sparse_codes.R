library(ProjectTemplate); load.project()
source("./munge/webster_depmap.R")

out_path <- file.path(".", "output", "06-sparse_codes")
create_output_folder(out_path)

# Functions ----------------------------------------------------------------



generate_sparse_code_heatmaps <- function(gene) { #cds_id expected
  tmp_gene_mat <- extract_atoms(webster_depmap %>% get_gene_mat(), gene)
  symbol <- gene %>% convert_cds_to_symbol()
  #Output 1: gene heatmap loadings
  
  pheatmap::pheatmap(tmp_gene_mat, main = gene, breaks = seq(-max(abs(tmp_gene_mat)),max(abs(tmp_gene_mat)), length.out =101),
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                     filename = file.path(out_path, sprintf("%s_gene_loadings.pdf",gene)), width = 6, height = 9)
  
  
  
  tmp_atom_mat <- webster_depmap %>% get_cell_mat() %>% set_colnames(paste("V", 1:webster_depmap$rank, sep = "")) %>%
    magrittr::extract(,colnames(tmp_gene_mat)[tmp_gene_mat[gene,] %>% order(decreasing = T)])
  

  #Output 2: fitness effects of functions
  atom_heatmap <- pheatmap::pheatmap(tmp_atom_mat, cluster_cols = F, show_rownames = F, show_colnames = T,
                                           color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(), 
                                           breaks = seq(-0.15, 0.15, length.out=101),
                                           cellwidth = 10,
                                           cellheight = 0.20, filename = file.path(out_path, sprintf("%s_dict_heatmap.pdf",gene)))
  
  
  #Output 3: fitness effects of genes, original and reconstructed. Preserve row order.
  tmp_dep <- cbind(fitness_mat_2[,gene], fitness_mat[,gene]) %>% 
    magrittr::extract(atom_heatmap$tree_row$order,)
  
  pheatmap::pheatmap(tmp_dep,
                     show_colnames = F, show_rownames = F, cluster_cols = F, cluster_rows = F,
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), 
                     breaks = seq(-4, 4, length.out=101),
                     cellwidth = 10,
                     cellheight = 0.20, filename = file.path(out_path, sprintf("%s_measured_recon_heatmap.pdf",gene)))
  
}

# Data --------------------------------------------------------------------


fitness_mat <- avana_19q4_webster %>% 
  t() %>% scale() %>% t()

fitness_mat_2 <- recon(webster_depmap) %>% t() %>% scale() %>% t()

# Heatmaps -------------------------------------------------------------------
#Used in paper
generate_sparse_code_heatmaps(convert_genes("SHOC2", 'symbol', 'cds_id'))

#Other examples
generate_sparse_code_heatmaps(convert_genes("C12orf49", 'symbol', 'cds_id'))
generate_sparse_code_heatmaps(convert_genes("C7orf26", 'symbol', 'cds_id'))


                   