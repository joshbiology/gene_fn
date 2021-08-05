library(ProjectTemplate); load.project()
library(igraph)


out_path <- "./output/02-genotoxic"

# Functions ---------------------------------------------------------------


generate_graph <- function(mat, rank = 5) {
  
  tmp <- cosine_sim(mat) %>% 
    edgeweight_symmetric_rank
  
  tmp[tmp >rank] <- NA
  
  tmp[upper.tri(tmp)] <- NA
  
  tmp[!is.na(tmp)] <- 1
  
  tmp2 <- igraph::graph_from_adjacency_matrix(tmp, mode = "undirected", weighted = NULL)
  
  return((tmp2))
  
}

graph_objective <- function(factorized_mat, orig_mat) {
  sum(diag(t(factorized_mat) %*% 
             (igraph::laplacian_matrix(generate_graph(orig_mat))) %*% 
             factorized_mat))
}

flip_max_cols <- function(source_mat) {
  flip_index <- abs(matrixStats::colMaxs(source_mat))-abs(matrixStats::colMins(source_mat))
  flip_index[flip_index < 0] <- -1
  flip_index[flip_index >= 0] <- 1
  
  #https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
  out <- source_mat %*% diag(flip_index) 
  rownames(out) <- rownames(source_mat)
  colnames(out) <- colnames(source_mat)
  return(out) 
}


#https://rdrr.io/cran/bigpca/src/R/bigpca.R
quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}


# Load genesets --------------------------------------------------------------------
essential_genes <- read_tsv("./data/raw/avana_v33-essential-genes.tsv")$gene %>% convert_genes("cds_id", "symbol")

tsg <- read_tsv("./data/raw/hart/tsg_table.txt")

olivieri_genes <- read_excel("./data/raw/durocher/cell/mmc4-2.xlsx", sheet = 1, skip = 1) %>%
  pivot_longer(names_to = "Pathway", values_to = "Gene", everything()) %>% 
  arrange(Pathway)

observed_clusters <- list(HIPPO = c("NF2", "FRYL", "AMOTL2", "LATS2", "TAOK1"),
  ESS = setdiff(essential_genes, olivieri_genes$Gene),
  TSG = tsg$GENE) %>% 
  enframe("Pathway", "Gene") %>% 
  unnest(Gene)

olivieri_genes <- olivieri_genes %>% rbind(observed_clusters) %>% na.omit()

chem_levels <- c("Cisplatin-3", "Cisplatin-2","Cisplatin-1", "UV", "Formaldehyde", "MNNG", "MMS", "CPT-2", "CPT-1", "Olaparib",  
                 "HU-acute", "HU-chronic", "MLN4924", "AZD6738", "Pyridostatin", "Doxorubicin", "Etoposide",   "IR",   "KBrO3",  "Bleomycin", "H2O2",
                 "CD437", "Gemcitabine", "BPDE", "illudinS", "Duocarmycin", "Calicheamicin", "Trabectedin", "PhenDC3", "ICRF-187")


# Load and filter fitness data --------------------------------------------


olivieri_mat <- read_excel("./data/raw/durocher/cell/mmc2-3.xlsx", sheet = 1, na = "NA") %>% 
  column_to_rownames("Gene") %>% 
  as.matrix


#Feature selection
olivieri_mat[is.na(olivieri_mat)] <- 0

stat_df <- tibble(Gene = rownames(olivieri_mat),
                  Var = matrixStats::rowVars(olivieri_mat),
                  Mean = matrixStats::rowMeans2(olivieri_mat),
                  Rank = order(order(desc(Var))))

stat_df %>% 
  ggplot(aes(Mean, Var)) +
  geom_point()

ggplot(stat_df, aes(sample = Var)) +
  geom_qq(size = 0.75) +
  geom_hline(yintercept = 3, color = "red") +
  ggsave(file.path(out_path, "genotoxic_gene_selection.png"), width = 4, height = 4)

damage_genes <- stat_df %>%
  filter(Var > 3) %>%
  pull(Gene)

# Exploring the fitness data --------------------------------------------
gene_annot <- stat_df %>% select(Gene, Var) %>% 
  left_join(olivieri_genes %>% mutate(Is_Member = 1) %>% unique()) %>% 
  pivot_wider(names_from = Pathway, values_from = Is_Member, values_fill = list(Is_Member = 0)) %>%
  select(-"NA") %>% 
  column_to_rownames("Gene")

umap_out <- umap::umap(olivieri_mat[damage_genes,], metric = "pearson", random_state = 2)

umap_out$layout %>% 
  as_tibble(rownames = "Gene") %>% 
  left_join(olivieri_genes) %>% 
  ggplot(aes(V1, V2, color = Pathway, label = Gene)) +
  geom_point() +
  ggsave(file.path(out_path, "umap_durocher_raw_data.pdf"), width = 5, height = 4)

library(tidygraph)
library(ggraph)

tidy_gene_graph <- as_tbl_graph(generate_graph(olivieri_mat[damage_genes,])) %>% 
  activate(nodes) %>% 
  left_join(olivieri_genes %>% group_by(Gene) %>% slice(1) %>% ungroup, by = c("name" = "Gene"))

ggraph(tidy_gene_graph) + 
  geom_edge_link(color = "gray20", alpha = 0.2) + 
  geom_node_point(size = 2.5, aes(colour = Pathway)) + 
  ggsave(file.path(out_path, "durocher_gene_graph.pdf"), width = 6, height = 5)

tidy_cl_graph <- as_tbl_graph(generate_graph(olivieri_mat[damage_genes,] %>% t()))

ggraph(tidy_cl_graph) + 
  geom_edge_link(color = "gray20", alpha = 0.2) + 
  geom_node_point(size = 2.5) + 
  theme_graph() +
  geom_node_text(aes(label = name), nudge_x = 0.1, nudge_y = -0.1) +
  ggsave(file.path(out_path, "durocher_cl_graph.pdf"), width = 6, height = 5)

# Run DGRDL ---------------------------------------------------------------

genotoxic_input <- olivieri_mat[damage_genes,] %>% t()

#For saving into the cache for use in other scripts:
#ProjectTemplate::cache("genotoxic_input")

export_for_matlab(olivieri_mat[damage_genes,] %>% t(), "./output/02-genotoxic/durocher_matlab.csv")

#Run DGRDL using written script

#Factorization output
mat_paths <- list("./data/interim/matlab/durocher", "./data/interim/matlab/durocher_no_graph") %>% 
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

grid_factorized <- map(mat_paths, function(x) {
  tmp <- R.matlab::readMat(x)
  tmp$X <- matrix(tmp$X, ncol = dim(tmp$X)[2])
  
  graphdl_to_factorized(import_graphdl(tmp), damage_genes, olivieri_mat %>% colnames())
})


params_df <- map_df(grid_factorized, function(x){
  tibble(K = attr(x, 'extras')$K %>% as.numeric,
         T_param = attr(x, 'extras')$T %>% as.numeric,
         Alpha = attr(x, 'extras')$alpha %>% as.numeric,
         Beta = attr(x, 'extras')$beta %>% as.numeric)
})

f_norms <- furrr::future_map_dbl(grid_factorized, function(x) {
  norm((olivieri_mat[damage_genes,] %>% t())-(x$cell_mat %*% t(x$gene_mat)), type = c("F"))
})

gene_laplacian <- furrr::future_map_dbl(grid_factorized, function(x) {
  graph_objective(x$gene_mat, olivieri_mat[damage_genes,])
})

cl_laplacian <- furrr::future_map_dbl(grid_factorized, function(x) {
  graph_objective(x$cell_mat, olivieri_mat[damage_genes,] %>% t())
})



# Process results ---------------------------------------------------------




metrics <- params_df %>% mutate(F_Norm = f_norms,
                                Gene_Laplacian = gene_laplacian,
                                Cell_Line_Laplacian = cl_laplacian
)  %>% 
  mutate(Graph_Regularized = case_when(Beta == 0 ~ F,
                                       T ~ T))

metrics %>% 
  ggplot(aes(K, F_Norm, color = factor(T_param), shape = Graph_Regularized)) +
  geom_smooth() +
  geom_point(size = 2) +
  labs(color='T') +
  facet_wrap(~Graph_Regularized, nrow = 1) +
  ggsave(file.path(out_path, "durocher_graphdl_fnorm.pdf"), width = 7, height = 3.5, device = cairo_pdf)

metrics %>% 
  ggplot(aes(K, Gene_Laplacian, color = factor(T_param), shape = Graph_Regularized)) +
  geom_point() +
  geom_smooth() +
  labs(color='T') +
  facet_wrap(~Graph_Regularized, nrow = 1) +
  scale_y_log10() +
  labs(y = "log(Gene Laplacian)") +
  ggsave(file.path(out_path, "durocher_graphdl_gene_lap.pdf"), width = 7, height = 3.5, device = cairo_pdf)

metrics %>% 
  ggplot(aes(K, Cell_Line_Laplacian, color = factor(T_param), shape = Graph_Regularized)) +
  geom_point() +
  geom_smooth() +
  labs(color='T') +
  facet_wrap(~Graph_Regularized, nrow = 1) +
  ggsave(file.path(out_path, "durocher_graphdl_cl_lap.pdf"), width = 7, height = 3.5, device = cairo_pdf)


# Set up Webster genotoxic  -----------------------------------------------------------------
source("./munge/webster_genotoxic.R")

# Portal export -----------------------------------------------------------


#Export for supplemental table
get_gene_mat(webster_genotoxic)[,function_order] %>% 
  set_colnames(function_names[function_order]) %>% 
  as_tibble(rownames = "Name") %>% 
  write_tsv(file.path(out_path, "durocher_inferred_functions.tsv"))

#Export for portal
webster_genotoxic %>% get_gene_mat() %>% 
  as_tibble(rownames = "Name") %>% 
  write_tsv(file.path(".", "output", "portal_assets", "durocher_gene_to_function.tsv"))

function_names %>% enframe("Name", "Display_Name") %>% mutate(Name = paste("V", Name, sep = "")) %>% 
  write_tsv(file.path(".", "output", "portal_assets", "durocher_fn_display_name.tsv"))


# Plots -------------------------------------------------------------------


#Fitness heatmaps
fitness_mat <- olivieri_mat[damage_genes,] %>% t()
fitness_mat <- scale(t(fitness_mat)) %>% t()



tmp <- pheatmap::pheatmap(fitness_mat, show_colnames = F, show_rownames = T, treeheight_row = 0, treeheight_col = 0,
                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), clustering_method = "ward.D2",
                          breaks = seq(-4, 4, length.out=101),
                          width = 7,
                          height = 4,
                          filename = file.path(out_path, "durocher_original_matrix_heatmap.pdf"))

pheatmap::pheatmap(fitness_mat, show_colnames = F, show_rownames = T, treeheight_row = 0, treeheight_col = 0,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), clustering_method = "ward.D2",
                   breaks = seq(-4, 4, length.out=101), annotation_col = gene_annot, annotation_legend = F,
                   width = 7,
                   height = 4 ,
                   filename = file.path(out_path, "durocher_original_matrix_heatmap_w_annot.pdf"))

#2021-4-1 output for figure generation
pheatmap::pheatmap(fitness_mat, show_colnames = T, show_rownames = T, treeheight_row = 0, treeheight_col = 0,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), clustering_method = "ward.D2",
                   breaks = seq(-3, 3, length.out=101),
                   cellwidth = 10,
                   cellheight = 10,
                   filename = file.path(out_path, "fixed_cell_width_durocher_original_matrix_heatmap.pdf"))

#Reconstructed

fitness_mat_2 <- recon(webster_genotoxic) 
fitness_mat_2 <- fitness_mat_2[tmp$tree_row$labels[tmp$tree_row$order], tmp$tree_col$labels[tmp$tree_col$order]]
fitness_mat_2 <- scale(t(fitness_mat_2)) %>% t()


rownames(fitness_mat_2) <- rownames(fitness_mat) 
colnames(fitness_mat_2) <- colnames(fitness_mat)

#Export reconstruction score
tibble(Name = damage_genes,
       Recon_Pearson = map_dbl(damage_genes, function(x) cor(olivieri_mat[x,], recon(webster_genotoxic)[,x]))) %>% 
  write_tsv(file.path(out_path, "durocher_recon_pearson.tsv"))

tmp2 <- pheatmap::pheatmap(fitness_mat_2,
                           show_colnames = F, show_rownames = T, cluster_cols = F, cluster_rows = F,
                           color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), 
                           breaks = seq(-4, 4, length.out=101),
                           width = 7,
                           height = 4,
                           filename = file.path(out_path, "durocher_reconstructed_matrix_heatmap.pdf"))

#2021-4-1 fixed width for figure generation
pheatmap::pheatmap(fitness_mat_2,
                   show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), 
                   breaks = seq(-3, 3, length.out=101),
                   cellwidth = 10,
                   cellheight = 10,
                   filename = file.path(out_path,"fixed_width_durocher_reconstructed_matrix_heatmap.pdf"))

#Gene factor heatmap
fitness_mat_3 <- webster_genotoxic$gene_mat %>% t() %>% set_colnames(damage_genes) %>% set_rownames(1:10) %>% 
  magrittr::extract(,tmp$tree_col$labels[tmp$tree_col$order])


pheatmap::pheatmap(fitness_mat_3[function_order,tmp$tree_col$labels[tmp$tree_col$order]], annotation_col = gene_annot,
                   show_colnames = F, show_rownames = T, cluster_cols = F, cluster_rows = F, annotation_legend = F,
                   color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), 
                   breaks = seq(-10, 10, length.out=101),
                   width = 7,
                   height = 4,
                   filename = file.path(out_path,"durocher_gene_loadings_heatmap.pdf"))

#Cell line heatmap
webster_genotoxic$cell_mat %>% set_rownames(colnames(olivieri_mat)) %>% set_colnames((1:10)) %>% 
  magrittr::extract(tmp$tree_row$labels[tmp$tree_row$order],function_order) %>% 
  pheatmap::pheatmap(cluster_cols = F, cluster_rows = F, show_colnames = F,
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(), 
                     breaks = seq(-0.75, 0.75, length.out=101),
                     width = 4,
                     height = 7,
                     filename = file.path(out_path,"durocher_cl_loadings_heatmap.pdf"))



# Export individual factor heatmaps ---------------------------------------


#Map one cell line factor for each 

map(1:10, function(x) {
  tmp_mat <- webster_genotoxic$cell_mat %>% set_rownames(colnames(olivieri_mat)) %>% set_colnames((1:10))
  
  
  tmp_mat[order((tmp_mat[,x])),x] %>% 
    matrix(ncol = 1) %>% 
    set_rownames(colnames(olivieri_mat)[order((tmp_mat[,x]))])  %>%
    pheatmap::pheatmap(cluster_cols = F, cluster_rows = F, show_rownames = T,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(), 
                       breaks = seq(-0.75, 0.75, length.out=101),
                       width = 1.85,
                       height = 5,
                       filename = sprintf(file.path(out_path,"durocher_cl_loadings_heatmap_%s.pdf"),function_names[x]))
  
})


#Map one gene factor for each 

map(1:10, function(x) {
  tmp_mat <- webster_genotoxic$gene_mat
  
  gene_selection <- order((tmp_mat[,x]), decreasing = T)[1:15]
  
  max_val <- max(tmp_mat[gene_selection,x]) + 10
  
  tmp_mat[gene_selection,x] %>% 
    matrix(ncol = 1) %>% 
    set_rownames(damage_genes[gene_selection])  %>%
    pheatmap::pheatmap(cluster_cols = F, cluster_rows = F, show_rownames = T, annotation_row = gene_annot, annotation_legend = F,
                       color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), 
                       breaks = seq(-max_val, max_val, length.out=101), main = function_names[x],
                       cellwidth = 10,
                       cellheight = 10,
                       filename = sprintf(file.path(out_path,"durocher_gene_loadings_heatmap%s.pdf"), function_names[x]))
  
})


#ATM ATR kinases

repair_kinases <- c("ATM", "ATR")

grid_factorized[[2]]$gene_mat %>% set_rownames(damage_genes) %>% magrittr::extract(repair_kinases,function_order) %>% 
  set_colnames(function_names[function_order]) %>% 
  pheatmap::pheatmap(cluster_cols = F, cluster_rows = F, show_rownames = T, show_colnames=T,
                     color = colorRampPalette(c("white", "#00AEEF"))(100), 
                     #breaks = seq(-0.75, 0.75, length.out=101),
                     width = 5,
                     height = 2.75,
                     filename = file.path(out_path,"durocher_atm_atr.pdf"))

# Joint UMAP embeding of fn and genes -----------------------------------------------------------

gene_df <- webster_genotoxic$gene_mat %>% t() %>% set_colnames(damage_genes) %>% set_rownames(1:10) %>% t() %>% as_tibble(rownames = "Gene") %>% 
  left_join(gene_annot %>% rownames_to_column("Gene"))

#Embed factors onto the UMAP
umap2_out <- rbind(olivieri_mat[damage_genes,],
                   webster_genotoxic$cell_mat %>% set_colnames(paste("V", 1:10, sep = "")) %>% t()) %>% 
  umap(metric = "pearson", n_neighbors = 15, random_state = 1)


umap2_df_gene <- umap2_out %>% 
  umap_to_df("Name") %>% 
  mutate(Type = c(rep("Gene", length(damage_genes)), rep("Factor", 10))) %>% 
  left_join(gene_df, by = c("Name"= "Gene"))

write_tsv(umap2_out %>%
            umap_to_df("Name") %>%
            mutate(Type = c(rep("Gene", length(damage_genes)), rep("Factor", 10))), file.path(out_path,"durocher_embedding_w_functions.tsv"))

umap2_df_gene %>% 
  ggplot(aes(V1, V2, color = Type, shape = Type)) +
  geom_point(alpha = 0.85) +
  scale_shape_manual(values=c(17, 16)) +
  scale_color_manual(values = c("Gene" = "gray90", "Factor" = "gray20")) +
  theme_void() +
  ggsave(file.path(out_path,"durocher_gene_fn_umap.pdf"), width = 2.75, height = 2.5, device = cairo_pdf)


umap2_df_gene %>% 
  pivot_longer(names_to = "Factor", values_to = "Loading", as.character(1:10)) %>% 
  mutate(Factor = factor(Factor, levels = function_order, labels = function_names[function_order])) %>% 
  ggplot(aes(V1, V2, color = Loading, shape = Type)) +
  geom_point(alpha = 0.85) +
  scale_shape_manual(values=c(17, 16)) +
  scale_color_gradient2(low = "#BD6C33", mid = "gray90", high =  "#00AEEF", midpoint = 0, limits=c(-10, 10), oob = scales::squish) +
  facet_wrap(~Factor) +
  theme_void() +
  ggsave(file.path(out_path,"durocher_gene_loadings_umap.pdf"), width = 8, height = 7, device = cairo_pdf)

umap2_df_gene %>% 
  left_join(olivieri_genes, by = c("Name" = "Gene")) %>% 
  ggplot(aes(V1, V2, color = Pathway, shape = Type)) +
  scale_shape_manual(values=c(17, 16)) +
  geom_point() +
  ggsave(file.path(out_path,"umap_durocher_with_lit_annotations.pdf"), width = 5, height = 4)


# Export Function annotations ---------------------------------------------
get_cell_mat(webster_genotoxic) %>% 
  set_colnames(function_names) %>% 
  t() %>% 
  as_tibble(rownames = "Name") %>% 
  write_tsv(file.path(out_path, "durocher_cell_mat.tsv"))

#For portal
get_cell_mat(webster_genotoxic) %>% 
  t() %>% 
  as_tibble(rownames = "Name") %>% 
  write_tsv(file.path(".","output", "portal_assets", "durocher_cell_mat.tsv"))


# functions for Ploting pleiotropy networks ----------------------------------------------------------------

threshold <- 0.25
loadings_df <- webster_genotoxic %>%
  get_gene_mat_stdev() %>% 
  set_colnames(function_names) %>% 
  as_tibble(rownames=  "Gene") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", -Gene) %>% 
  filter(abs(Loading) > threshold)



factor_pleiotrop_df <- loadings_df %>% 
  left_join(loadings_df, by = "Gene") %>% 
  filter(Function.x != Function.y, Function.x < Function.y) %>% 
  left_join(tibble(Gene = damage_genes,
                   Recon_Pearson = map_dbl(damage_genes, function(x) cor(olivieri_mat[x,], recon(webster_genotoxic)[,x])))) %>% 
  filter(Recon_Pearson > 0.4)


nrow(factor_pleiotrop_df)

plot_pleiotropy_network <- function(fns, pathway_name) {
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




# Interpolation -----------------------------------------------------------
mock_rad51b <- genotoxic_input[,"H2AFX"]/norm(genotoxic_input[,"H2AFX"], type = "2") - get_cell_mat(webster_genotoxic)[,7] + get_cell_mat(webster_genotoxic)[,8]
cor(mock_rad51b, genotoxic_input[,"RAD51B"])

qplot(genotoxic_input[,"RAD51B"], genotoxic_input[,"H2AFX"], ) +
  labs(y = "Measured H2AFX fitness", x = "Measured RAD51B fitness") + 
  ggsave(file.path(out_path, "rad51_measured.pdf"), width = 4, height = 4, device = cairo_pdf)

qplot(genotoxic_input[,"RAD51B"], mock_rad51b) +
  labs( x = "Measured RAD51B fitness", y = "Measured H2AFX - End Joining + Fanconi Anemia")  + 
  ggsave(file.path(out_path, "rad51_artificial.pdf"), width = 4, height = 4, device = cairo_pdf)

#heatmaps


# Run other factorizations ------------------------------------------------

#PCA
pc <- prcomp(genotoxic_input %>% t())
pc_dim <- quick.elbow(pc$sdev^2)

plot(pc, type = "lines", npcs = 31)

pca_factor <- pca_to_factorized(pc, pc_dim)


#ICA
K_param<- seq(2, 30, 1)

#Make genes the S matrix so that the non-Gaussianiity is imposed on the profiles over cell line measurements
ics <- purrr::map(K_param, ~icasso_ica(olivieri_mat[damage_genes,],nbComp=., alg.type="deflation", nbIt=20,
                                       funClus="hclust", method="average"))

mstd_df <- purrr::map(ics, ~.$Iq %>% as.data.frame() %>% 
                        set_colnames("Stability") %>% 
                        mutate(IC = row_number())) %>% 
  set_names(K_param) %>% 
  enframe("Rank") %>%
  mutate(Rank = as.numeric(Rank)) %>% 
  unnest(value) 

mstd_df %>% 
  ggplot(aes(IC, Stability, group = Rank, color = Rank)) +
  geom_line() +
  geom_point() +
  ylim(c(0, 1)) +
  scale_color_viridis_c() +
  ggtitle("MSTD plot for Durocher ICA") +
  geom_vline(xintercept = 13, color = "red", linetype = "dashed") +
  ggsave(file.path(out_path, "mstd_durocher.pdf"), width = 5, height = 4)

ica_factor <- fastICA_to_factorized(ics[[which(K_param == 10)]]) 



# Perform geneset operations ----------------------------------------------

one_hot_mat <- olivieri_genes %>% 
  filter(Pathway != "HIPPO", Pathway != "OTHER") %>% 
  mutate(Is_Pathway = T) %>% 
  pivot_wider(names_from = "Pathway", values_from = "Is_Pathway", values_fill = F) %>% 
  column_to_rownames("Gene") %>% 
  as.matrix() %>% 
  "*"(1) #trick to convert logical to numeric


plot_roc_over_loadings <- function(factor_obj) {
  gene_loadings <- factor_obj %>% 
    get_gene_mat() %>% 
    as_tibble(rownames = "Gene") %>% 
    pivot_longer(names_to = "Factor", values_to = "Loading", -Gene)
  
  out <- map_dfr(1:ncol(one_hot_mat), function(x) {
    gene_loadings %>% 
      inner_join(one_hot_mat[,x] %>% enframe("Gene", "Truth")) %>% 
      mutate(Truth = factor(Truth)) %>% 
      group_by(Factor) %>% 
      yardstick::roc_auc(Truth, Loading, event_level = "second")}, .id = "Index")   #https://yardstick.tidymodels.org/

  out %>% 
    left_join(tibble(Index = 1:ncol(one_hot_mat) %>% as.character, Name = colnames(one_hot_mat))) %>% 
    ggplot(aes(x = Name, y = .estimate, color = Factor)) +
    geom_hline(yintercept = 0.5, linetype='dotted') +
    geom_point() +
    ylim(c(0, 1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  
    
}

plot_roc_over_loadings(webster_genotoxic) +
  ggsave(file.path(out_path, "roc_geneset_webster.pdf"), width = 4, height = 4, device = cairo_pdf)

plot_roc_over_loadings(pca_factor) +
  ggsave(file.path(out_path, "roc_geneset_pca.pdf"), width = 4, height = 4, device = cairo_pdf)

plot_roc_over_loadings(ica_factor) +
  ggsave(file.path(out_path, "roc_geneset_ica.pdf"), width = 4, height = 4, device = cairo_pdf)


