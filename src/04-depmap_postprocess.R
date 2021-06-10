
library(ProjectTemplate); load.project()

require(future)
require(purrr)
plan(multisession)


out_path <- file.path(".", "output", "04-depmap_postprocess")
# Functions ---------------------------------------------------------------



graph_objective_precomp <- function(factorized_mat, L_mat) {
  sum(diag(Rfast::Crossprod(factorized_mat,L_mat) %>% 
             Rfast::mat.mult(factorized_mat)))
}

generate_graph <- function(mat, rank = 5) {
  require(igraph)
  tmp <- cosine_sim(mat) %>% 
    edgeweight_symmetric_rank
  
  tmp[tmp >rank] <- NA
  
  tmp[upper.tri(tmp)] <- NA
  
  tmp[!is.na(tmp)] <- 1
  
  tmp2 <- igraph::graph_from_adjacency_matrix(tmp, mode = "undirected", weighted = NULL)
  
  return((tmp2))
  
}


import_graphDL  <- function(x) {
  tmp <- R.matlab::readMat(x)

  graphdl_to_factorized(import_graphdl(tmp), colnames(avana_19q4_webster), rownames(avana_19q4_webster))
}

extract_params <- function(output) {
  map_dfr(output, function(x){
    tibble(K = attr(x, 'extras')$K %>% as.numeric,
           T_param = attr(x, 'extras')$T %>% as.numeric,
           Alpha = attr(x, 'extras')$alpha %>% as.numeric,
           Beta = attr(x, 'extras')$beta %>% as.numeric,
           Seed = attr(x, 'extras')$seed %>% as.numeric,
           Num_Neigh_Gene = attr(x, 'extras')$num.neighbor.gene[1])
  })
  
}

generate_metrics <- function(output) {
  params_df <- extract_params(output)
  
  f_norms <- furrr::future_map_dbl(output, function(x) {
    norm(avana_19q4_webster-recon(x), type = c("F"))
  })
  
  
  g_laps <- furrr::future_map_dbl(output, function(x) {
    graph_objective_precomp(x$gene_mat, gene_laplacian)
  })
  
  params_df %>% 
    mutate(F_Norm = f_norms,
           Gene_Lap  =g_laps)
}

# Data --------------------------------------------------------------------


load("./cache/avana_19q4_webster.RData")
gene_laplacian <- igraph::laplacian_matrix(generate_graph(avana_19q4_webster%>% t())) %>% as.matrix()


# Grid --------------------------------------------------------------------
#Grid factorization output: many values of K, many values of T, one random seed
mat_paths <- list("./data/interim/matlab/depmap_grid") %>% 
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

depmap_grid_output <- map(mat_paths, import_graphDL)


tictoc::tic()
depmap_grid_metrics_df <- generate_metrics(depmap_grid_output)
tictoc::toc()

write_tsv(depmap_grid_metrics_df, path = file.path(out_path, "depmap_grid_metrics.tsv"))
#If loading from saved TSV:
#depmap_grid_metrics_df <- read_tsv(file.path(out_path, "depmap_grid_metrics.tsv"))

#Save space
rm(depmap_grid_output)

depmap_grid_metrics_df %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>% 
  filter(K < 600) %>% 
  
  ggplot(aes(x = K, y = Value, group = T_param, color = factor(T_param))) +
  geom_point() +
  geom_line() +
  facet_wrap(~Metric, scales ="free_y", ncol = 2) +
  ggsave(file.path(out_path,"depmap_grid_metrics.pdf"), width = 10, height = 6, device = cairo_pdf)


# Wide ------------------------------------------------------------------

#Wide factorization output: many values of K, one value of T, one random seed
mat_paths <- list("./data/interim/matlab/depmap_wide") %>% 
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

depmap_wide_output <- map(mat_paths, import_graphDL)

tictoc::tic()
depmap_wide_metrics_df <- generate_metrics(depmap_wide_output)
tictoc::toc()

write_tsv(depmap_wide_metrics_df, path = file.path(out_path, "depmap_wide_metrics.tsv"))

depmap_wide_metrics_df %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>% 
  ggplot(aes(x = K, y = Value)) +
  geom_jitter() +
  geom_smooth() +
  geom_vline(xintercept = 220) +
  xlim(c(0, 600)) +
  facet_grid( Metric ~Num_Neigh_Gene, scales ="free_y") +
  ggsave(file.path(out_path,"depmap_wide_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)

depmap_wide_metrics_df %>% 
  arrange(K) %>% 
  mutate(Diff_F = c(0,diff(F_Norm)),
         Diff_Gene_L = c(0,diff(Gene_Lap))) %>% 
  slice(-1) %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", c("Diff_F", "Diff_Gene_L")) %>% 
  filter(K < 500) %>% 
  ggplot(aes(x = K, y = Value))+
  geom_vline(xintercept = 220) +
  xlim(c(0, 500)) +
  geom_jitter() +
  geom_smooth() +
  facet_grid( Metric ~Num_Neigh_Gene, scales ="free_y") +
  ggsave(file.path(out_path,"depmap_wide_marginal_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)

#Save space
rm(depmap_wide_output)


# Deeper dive -------------------------------------------------------------


#Deep factorization output: few values of K, one value of T, many random seeds
mat_paths <- list("./data/interim/matlab/depmap_deep") %>% 
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

depmap_deep_output <- map(mat_paths, import_graphDL)

tictoc::tic()
depmap_deep_metrics_df <- generate_metrics(depmap_deep_output)
tictoc::toc()

write_tsv(depmap_deep_metrics_df, path = file.path(out_path, "depmap_deep_metrics.tsv"))

depmap_deep_metrics_df %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>% 
  ggplot(aes(x = K, y = Value)) +
  geom_jitter() +
  geom_smooth() +
  facet_grid( Metric ~Num_Neigh_Gene, scales ="free_y") +
  ggsave(file.path(out_path,"depmap_deep_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)

grouped_metrics <- depmap_deep_metrics_df %>% 
  group_by(K) %>% summarize(Mean_F_Norm = mean(F_Norm), Mean_Gene_Lap = mean(Gene_Lap))
  
grouped_metrics %>% 
  mutate(Diff_F = c(0,diff(Mean_F_Norm)),
         Diff_Gene_L = c(0,diff(Mean_Gene_Lap))) %>% 
  slice(-1) %>% 
  pivot_longer(names_to = "Metric", values_to = "Value", c("Diff_F", "Diff_Gene_L")) %>% 
  filter(K < 500) %>% 
  ggplot(aes(x = K, y = Value))+
  geom_jitter() +
  geom_smooth() +
  facet_wrap(  ~Metric, ncol = 1, scales ="free_y") +
  ggsave(file.path(out_path,"depmap_deep_marginal_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)



# Inidividual gene examination -------------------------------------------------------
#returns the genes and columns invovled.
extract_atoms <- function(mat, gene, loading_threshold = 8) {
  # set col names
  colnames(mat) <- paste("V", 1:ncol(mat), sep = "")
  
  gene_loadings <- mat[gene, ]
  col_index <- abs(gene_loadings) > 0
  
  row_tmp <- rowSums(abs(mat[,col_index]))
  row_index <- row_tmp > loading_threshold
  
  
  return(mat[row_index, col_index])
}


generate_seed_heatmaps <- function(gene, output= depmap_deep_output) {
  
  require(RColorBrewer)
  params_df <- output %>% extract_params()
  heatmap_names <- map(1:length(output), function(x) sprintf("%s_nneigh=%d_k=%d_seed=%d_T=%d_Alpha=%f", 
                                                              gene,
                                                              params_df$Num_Neigh_Gene[x],
                                                              params_df$K[x],
                                                              params_df$Seed[x],
                                                              params_df$T_param[x],
                                                              params_df$Alpha[x]))
  

  
  
  heatmaps <- walk(1:length(output), function(x) {
    tmp_mat <- output[[x]]$gene_mat %>% set_rownames(attr( output[[1]], 'gene_names'))
    atoms <- extract_atoms(tmp_mat, gene)
    pheatmap::pheatmap(atoms, main = heatmap_names[x], breaks = seq(-max(abs(atoms)),max(abs(atoms)), length.out =101),
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                       filename = sprintf("./output/04-depmap_postprocess/%s.pdf",heatmap_names[x]), width = 6, height = 9)
  })
  
}


genes = c(shoc2 = "SHOC2",
          kras = "KRAS",
          c7orf26 = "C7orf26",
          ints6 = "INTS6",
          ints12 = "INTS12",
          tp53 = "TP53", 
          mtor = "MTOR",
          braf = "BRAF",
          hmgcr = "HMGCR",
          mvk = "MVK",
          ubiad = "UBIAD1",
          smarca4= "SMARCA4",
          nrf2 = "NFE2L2",
          cdk6 = "CDK6",
          skp2 = "SKP2",
          egfr = "EGFR",
          rictor = "RICTOR",
          tgfb = "TGFBR1",
          ctnnb1 = "CTNNB1",
          hippo = "NF2",
          myc = "MYC",
          max = "MAX",
          emsy = "EMSY",
          ncor = "NCOR1",
          wrn = "WRN",
          kat2a = "KAT2A")



other_genes = c("TXNRD1", "PRKRA", "RPP25L", "WDR1",  "LEMD3", "EGLN1", "VRK1", "BIRC6", "ACSL3")#"TMPO",

convert_genes(genes, 'symbol', 'cds_id')


convert_genes(genes, 'symbol', 'cds_id') %in% (avana_19q4_webster %>% colnames())
convert_genes(other_genes, 'symbol', 'cds_id') %in% (avana_19q4_webster %>% colnames())


#Goal here is to verify the factors used in the paper
walk(convert_genes(genes, 'symbol', 'cds_id'), generate_seed_heatmaps)
walk(convert_genes(other_genes, 'symbol', 'cds_id'), generate_seed_heatmaps)

generate_seed_heatmaps(convert_genes("GPX4",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("C12orf49",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("ACSL3",'symbol', 'cds_id'))

generate_seed_heatmaps(convert_genes("PEX26",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("SOS1",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("BRD9",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("BRD7",'symbol', 'cds_id'))

generate_seed_heatmaps(convert_genes("SMARCB1",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("ALG5",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("DHODH",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("NAMPT",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("DTYMK",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("DTYMK",'symbol', 'cds_id'))


generate_seed_heatmaps(convert_genes("KAT2A",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("TADA2A",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("TADA2B",'symbol', 'cds_id'))

generate_seed_heatmaps(convert_genes("VPS29",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("PRKRA",'symbol', 'cds_id'))
generate_seed_heatmaps(convert_genes("EP300",'symbol', 'cds_id'))


