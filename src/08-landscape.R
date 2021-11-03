
#Visualizing graphDL output
library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "08-landscape")
create_output_folder(out_path)


# Functions ---------------------------------------------------------------

geneset_specificity <- function(factor_genes, geneset, all_genes) {

  return(mltools::mcc(preds = all_genes %in% factor_genes,
                      actuals = all_genes %in% geneset))

}


gene_mat_to_set <- function(gene_mat, gene_names) {
  #Takes any non-negative loading for a function and creates a geneset
  map(1:ncol(gene_mat), ~gene_names[gene_mat[,.] != 0])
}


softmax <- function(par){
  #https://rpubs.com/FJRubio/softmax
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk)))
  }
  val <- exp(par - Lk)
  return(val)
}


# Gene and Factor UMAP ----------------------------------------------------

source("./munge/webster_depmap.R")

combo <- cbind(webster_depmap %>% get_cell_mat %>% set_colnames(paste("V", 1:ncol(.), sep = "")), avana_19q4_webster)

depmap_umap <- umap(combo %>% t(), metric = 'pearson', n_neighbors = 10, random_state =1.1)

ProjectTemplate::cache("depmap_umap")

plotting_df <- depmap_umap$layout %>%
  as_tibble(rownames = "Name") %>%
  mutate(Type = c(rep("Function", ncol(get_cell_mat(webster_depmap))), rep("Gene", dim(avana_19q4_webster)[2]))) %>%
  mutate(Gene_Name = convert_genes(Name, from = "cds_id", "symbol")) %>%
  rename(X = V1, Y = V2) %>%
  mutate(Y = -Y)

ProjectTemplate::cache("plotting_df")

#For reproducibility:
#load("./cache/plotting_df.RData")

write_tsv(plotting_df, path = file.path(out_path, "depmap_embedding.tsv"))

# Individual plots --------------------------------------------------------


#Joint embedding
plotting_df %>%
  ggplot(aes(X, Y)) +
  geom_point(data = subset(plotting_df, Type == "Gene"), size = 0.5, alpha = 0.3, color = "gray70") +
  geom_point(data = subset(plotting_df, Type == "Function"),size = 1, shape = 17)
  ggsave(filename = file.path(out_path, "depmap_embedding.pdf"), width = 5, height = 5, device = cairo_pdf)

#Joint embedding void teme
plotting_df %>%
  ggplot(aes(X, Y)) +
  geom_point(data = subset(plotting_df, Type == "Gene"), size = 0.5, alpha = 0.3, color = "gray80") +
  geom_point(data = subset(plotting_df, Type == "Function"),size = 1, shape = 17) +
  theme_void()
  ggsave(filename = file.path(out_path, "depmap_embedding_void.pdf"), width = 4, height = 4, device = cairo_pdf)

#Genes only
plotting_df %>%
  ggplot(aes(X, Y)) +
  geom_point(data = subset(plotting_df, Type == "Gene"), size = 0.5, alpha = 0.3, color = "gray70")
  ggsave(filename = file.path(out_path, "depmap_embedding_gene.pdf"), width = 5, height = 5, device = cairo_pdf)

#Fns only
plotting_df %>%
  ggplot(aes(X, Y)) +
  geom_point(data = subset(plotting_df, Type == "Function"),size = 1, shape = 17)
  ggsave(filename = file.path(out_path, "depmap_embedding_fn.pdf"), width = 5, height = 5, device = cairo_pdf)

#Fn labeled
plotting_df %>%
  ggplot(aes(X, Y)) +
  geom_point(data = subset(plotting_df, Type == "Function"),size = 1, shape = 17) +
  geom_text(data = subset(plotting_df, Type == "Function"), aes(label = Name), size = 2, nudge_y = 0.15)
  ggsave(filename = file.path(out_path, "depmap_embedding_fn_labeled.pdf"), width = 5, height = 5, device = cairo_pdf)

#Genes by median essentiality
plotting_df %>%
  left_join(tibble(Name = colnames(avana_19q4_webster),
                   Median = avana_19q4_webster %>% matrixStats::colMedians())) %>%
  filter(Type == "Gene") %>%
  ggplot(aes(X, Y)) +
  geom_point(size = 0.5, alpha = 0.75, aes(color = Median)) +
  scale_color_distiller(palette = "RdBu")+
  theme(legend.position = "bottom")
  ggsave(filename = file.path(out_path, "depmap_embedding_gene_median_essential.pdf"), width = 4, height = 4.5, device = cairo_pdf)

#Genes by Webster reconstruction
plotting_df %>%
  left_join(gene_recon_webster, by = c("Name" = "Gene")) %>%
  filter(Type == "Gene") %>%
  ggplot(aes(X, Y)) +
  geom_point(size = 0.5, alpha = 0.3, aes(color = Recon_Pearson)) +
  scale_color_viridis_c() +
  labs(color = "Pearson") +
  theme(legend.position = "bottom")
  ggsave(filename = file.path(out_path, "depmap_embedding_gene_recon.pdf"), width = 4, height = 4.5, device = cairo_pdf)

#Genes by original variance
plotting_df %>%
  left_join(tibble(Name = colnames(avana_19q4_webster),
                   Var = avana_19q4_webster %>% matrixStats::colVars())) %>%
  filter(Type == "Gene") %>%
  ggplot(aes(X, Y)) +
  geom_point(size = 0.5, alpha = 0.3, aes(color = log(Var))) +
  scale_color_viridis_c(option="magma") +
  theme(legend.position = "bottom")
  ggsave(filename = file.path(out_path, "depmap_embedding_gene_var.pdf"), width = 4, height = 4.5, device = cairo_pdf)

# Project Localization data into Webster --------------------------------------------------------------------

source("./munge/subcell.R")

loc_genes <- intersect(rownames(get_gene_mat(webster_depmap)), bioid_loc$entrezgene %>% convert_genes("entrez_id", "cds_id"))

#NMF localization profiles
nmf_projected <- crossprod(get_gene_mat(webster_depmap)[loc_genes,] %>% abs, bioid_nmf[loc_genes %>% convert_cds_to_symbol,]) %>%
  aaply(1, function(x)x/sum(x)) %>%
  set_rownames(paste("V", 1:webster_depmap$rank, sep = ""))

entropy_localization <-nmf_projected %>% aaply(1, entropy::entropy)

nmf_cluster <- pheatmap::pheatmap(nmf_projected,
                   cluster_cols = T, cluster_rows = T, clustering_distance_cols = "correlation", labels_col = bioid_meta$Location %>% str_sub(1,25),
                   color = viridis::inferno(100), annotation_row = entropy_localization %>% enframe(value = "Entropy") %>% column_to_rownames("name"),
                   filename = file.path(out_path, "function_localization_profile.pdf"), width = 10, height = 9)


ProjectTemplate::cache('nmf_projected')

# Curate localization data -------------------------------------------------

max_nmf_df <- nmf_projected %>%
  adply(1, function(x) {
    index = which(max(x) == x)
    tibble(NMF_Column_Name = colnames(nmf_projected)[index], Score = x[index])
    }) %>%
  rename("Name" = "X1") %>%
  as_tibble() %>%
  left_join(entropy_localization %>% enframe("Name", "Entropy")) %>%
  left_join(bioid_meta) %>%
  mutate(Specificity = scales::rescale(Entropy, to=c(1,0))) %>%
  mutate(Location = factor(Location, levels= bioid_meta$Location[nmf_cluster$tree_col$order]))

ProjectTemplate::cache("max_nmf_df")

# Graph dendrogram -------------------------------------------------------


bioid_colors <- tibble(rank =  nmf_cluster$tree_col$order,
                    Color = c(rep("#21409a", 4),
                              rep("#ed1c24", 3),
                              rep("#2f4f4f", 2),
                              rep("#ff8c00", 2),
                              rep("#009444", 4),
                              rep("#ff69b4", 2),
                                rep("#00aeef", 3)),
                    Compartment = c(rep("Nucleus", 4),
                                    rep("Cytoskeleton", 3),
                                    rep("Misc.", 2),
                                    rep("Mito", 2),
                                    rep("ER", 4),
                                    rep("Membrane", 2),
                                    rep("Trafficking", 3)))

ProjectTemplate::cache("bioid_colors")


require(ggraph)
require(tidygraph)

loc_graph <- as.dendrogram(nmf_cluster$tree_col) %>% as_tbl_graph %>%
  left_join(bioid_meta, by = c("label" = "NMF_Column_Name")) %>%
  left_join(bioid_colors)

loc_graph %>%
  ggraph('dendrogram') +
  geom_edge_diagonal() +
  geom_node_text(aes(label = Location, color = Color), angle = 90, hjust=1, nudge_y = -0.15) +
  geom_node_point(aes(color = Color)) +
  scale_colour_identity() +
  ylim(-10, NA)
  ggsave(file.path(out_path, "loc_dendrogram.pdf"), width = 5, height = 8, device = cairo_pdf)

# Plot embedding ----------------------------------------------------------
#Plotting locations
loc_plot_df <- plotting_df %>%
  left_join(max_nmf_df) %>%
  left_join(bioid_colors)


loc_plot <- loc_plot_df %>%
  filter(Type == "Function") %>%
  ggplot(aes(X, Y)) +
  geom_point(data = loc_plot_df %>% filter(Type == "Gene"), size = 0.5, alpha = 0.4, color = "gray90") +
  geom_point(shape = 17, aes(color = Color, alpha = Specificity)) +
  scale_colour_identity() +
  scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
  theme_void()

loc_plot +
  ggsave(file.path(out_path, "loc_embedding.pdf"), width = 6, height = 5)


loc_plot_df %>%
  filter(Type == "Function") %>%
  ggplot(aes(X, Y)) +
  geom_point(shape = 17, size= 0.75, aes(color = Color, alpha = Specificity)) +
  scale_colour_identity() +
  scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
  facet_wrap(~Color, nrow = 1) +
  theme(legend.position = "none")
  ggsave(file.path(out_path, "facet_loc_embedding.pdf"), width = 8, height = 1.75)



loc_plot_df %>%
  filter(Type == "Function") %>%
  ggplot(aes(X, Y)) +
  geom_point(data = loc_plot_df %>% filter(Type == "Gene"), size = 0.5, alpha = 0.4, color = "gray90") +
  geom_point(shape = 17, aes(color = Color, alpha = Specificity)) +
  geom_text(aes(label = Name, color = Color, alpha = Specificity)) +
  scale_colour_identity() +
  scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
  theme_void()
  ggsave(file.path(out_path, "loc_embedding_w_names.pdf"), width = 7, height = 6)

g1 <- loc_plot_df %>%
  ggplot(aes(X, Y, label = Name)) +
  geom_point(data = loc_plot_df %>% filter(Type == "Gene"), size = 0.5, alpha = 0.4, color = "gray50") +
  geom_point(data = loc_plot_df %>% filter(Type == "Function"), shape = 17, aes(color = Color, alpha = Specificity)) +
  scale_colour_identity() +
  scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
  theme_void()

plotly::ggplotly(g1)


# Plotting embedding without genes ---------------------------------------


fns_only <- umap::umap(webster_depmap %>% get_cell_mat %>% t(),  metric = 'pearson', n_neighbors = 10, random_state =1.1)

g3 <- fns_only %>%
  umap_to_df("Name") %>%
  left_join(max_nmf_df) %>%
  left_join(bioid_colors) %>%
  ggplot(aes(V1, V2, color = Color, alpha = Specificity, label=Name)) +
  geom_point(shape = 17) +
  scale_color_identity() +
  scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +theme_void()


g3 +
  ggsave(file.path(out_path, "fn_only.pdf"), width = 3, height = 2)
plotly::ggplotly(g3)

# Functions for Plotting frames ---------------------------------------------------------


plotting_window <- function(fn_names, nudge = 0.3) {

  tmp <- loc_plot_df %>% filter(Name %in% fn_names)
  window_x <- c(min(tmp$X)-nudge, max(tmp$X)+nudge)
  window_y <- c(min(tmp$Y)-nudge, max(tmp$Y)+nudge)
  loc_plot_df %>%
    filter(X > window_x[1], X < window_x[2], Y > window_y[1], Y < window_y[2])
}

plot_loadings_panel <- function(df, focus_fn, gene_labs = T){


  gene_loadings_df <- webster_depmap %>%
    get_gene_mat() %>%
    as_tibble(rownames = "Name")   %>%
    pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V")) %>%
    filter(Function == focus_fn)

  df <-   df %>%
    left_join(gene_loadings_df)

  if(gene_labs) {
    df %>%
      ggplot(aes(X, Y, color = Loading, shape = Type)) +
      geom_point(data = subset(df, Type == "Gene" & Loading == 0),alpha = 0.85) +
      geom_point(data = subset(df, Type == "Gene" & Loading != 0),alpha = 0.85) +
      geom_text(data = df %>% filter(Type == "Gene"), aes(label = Gene_Name), alpha = 0.5, color = "black", position = position_nudge(y = -0.01)) +

      geom_point(data = subset(df, Type == "Function"),alpha = 0.85) +
      geom_text(data = df %>% filter(Type == "Function"), aes(label = Name),alpha = 0.5, color = "black",position = position_nudge(y = -0.01)) +

      scale_shape_manual(values=c(17, 16)) +
      scale_color_gradient2(low = "#BD6C33", mid = "gray90", high =  "#00AEEF", midpoint = 0, limits=c(-10, 10), oob = scales::squish) +
      theme_void() +
      theme(legend.position = "NA")
  }
  else
    df %>%
    ggplot(aes(X, Y, color = Loading, shape = Type)) +
    geom_point(data = subset(df, Type == "Gene" & Loading == 0),alpha = 0.2, color = "gray90") +
    geom_point(data = subset(df, Type == "Gene" & Loading != 0),alpha = 0.85) +

    geom_point(data = subset(df, Type == "Function" & Name == focus_fn)) +

    scale_shape_manual(values=c(17, 16)) +
    scale_color_gradient2(low = "#BD6C33", mid = "gray90", high =  "#00AEEF", midpoint = 0, limits=c(-10, 10), oob = scales::squish) +
    theme_void() +
    theme(legend.position = "NA")

}


plot_gene_panel <- function(df){

  df %>%
    ggplot(aes(X, Y, label = Gene_Name)) +
    geom_point(data = df %>% filter(Type == "Gene"), size = 1.5, alpha = 0.5, color = "gray70") +
    geom_text(data = df %>% filter(Type == "Gene"), alpha = 0.5, color = "black") +
    geom_point(data = df %>% filter(Type == "Function"), size = 3,shape = 17, aes(color = Color, alpha = Specificity)) +
    scale_colour_identity() +
    scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
    theme_void()
}

plot_fn_panel <- function(df) {
  df %>%
    ggplot(aes(X, Y, label = Name)) +
    geom_point(data = df %>% filter(Type == "Gene"), size = 1.5, alpha = 0.5, color = "gray70") +
    geom_point(data = df %>% filter(Type == "Function"), size = 3,shape = 17, aes(color = Color, alpha = Specificity)) +
    geom_text(data = df %>% filter(Type == "Function"),  aes(color = Color, alpha = Specificity), nudge_y = .05) +
    scale_colour_identity() +
    scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
    theme_void()
}


# functions for Ploting pleiotropy networks ----------------------------------------------------------------


loadings_df <- webster_depmap %>%
  get_gene_mat() %>%
  as_tibble(rownames=  "Gene") %>%
  pivot_longer(names_to = "Function", values_to = "Loading", -Gene) %>%
  filter(abs(Loading) > 0)

threshold <- 0.05

factor_pleiotrop_df <- loadings_df %>%
  left_join(loadings_df, by = "Gene") %>%
  filter(Function.x != Function.y, Function.x < Function.y) %>%
  mutate(Fraction = abs(Loading.x)/(abs(Loading.x)+abs(Loading.y))) %>%
  filter(Fraction < 1-threshold, Fraction > threshold, Loading.x >1, Loading.y > 1) %>%
  arrange(Fraction)

plot_pleiotropy_network <- function(fns, pathway_name) {
  require(igraph)
  require(tidygraph)
  require(ggraph)
  tmp_graph <- as_tbl_graph( igraph::graph_from_data_frame(factor_pleiotrop_df %>%
                                                             filter(Function.x %in% fns & Function.y%in% fns) %>%
                                                             count(Function.x, Function.y) %>%
                                                             arrange(desc(n)), directed = F))

  ggraph(tmp_graph, layout = "circle") +
    geom_edge_link(color = "gray90", alpha = 1, aes(edge_width = n)) +
    geom_node_text(aes(label = name), size = 5, color = "black") +
    geom_node_point(size = 2.5, color = "black")
  #ggsave(paste("./output/figures/", pathway_name, "pleiotropy.pdf", sep = ""), width = 4.5, height = 2.5, device = cairo_pdf)

}




# Focus on Sterol / Peroxisome --------------------------------------------
sterol <- "V30"
perox <- "V51"


sterol_window <- plotting_window(c(sterol, perox), nudge = 0.05)

sterol_window %>%
  plot_gene_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "sterol_window.pdf"), width = 3, height = 3)

plot_pleiotropy_network(c(sterol, perox))


# Focus on Trafficking  --------------------------------------------
traf_expand_window <- plotting_window(c(sterol, perox), nudge = 1)

traf_expand_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "sterol_expand.pdf"), width = 3, height = 3)


# Focus on Retrograde --------------------------------------------

wash <- "V11"
hops <- "V169"
lamtor <- "V17"

retro_window <- plotting_window(c(wash, hops, lamtor), nudge =0.2)
retro_window  %>%
  plot_gene_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "retro_window.pdf"), width = 3, height = 3)

plot_pleiotropy_network(c(wash, hops, lamtor))


# Focus on Recycling --------------------------------------------
retro_expand_window <- plotting_window(c(wash, hops, lamtor), nudge =0.75)
retro_expand_window  %>%
  plot_fn_panel()  +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "retro_expand.pdf"), width = 3, height = 3)


# Focus on SCAR/WAVE  --------------------------------------------
scar_wave <- "V52"
focal <- "V6"

scar_window <- plotting_window(c(scar_wave, focal), nudge = 0.05)

scar_window%>%
  plot_gene_panel()  +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "scar_window.pdf"), width = 3, height = 3)

plot_pleiotropy_network(c(scar_wave, focal))


# Focus on membrane -------------------------------------------------------

scar_expanded_window <- plotting_window(c(scar_wave, focal), nudge = 0.5)

scar_expanded_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "scar_expanded.pdf"), width = 3, height = 3)



# Mito window  -------------------------------------------------------

mito_window <- plotting_window(c("V2", "V145"), nudge = 0.05)

mito_window %>%
  plot_gene_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "mito_window.pdf"), width = 3, height = 3)

# Mito expanded  -------------------------------------------------------

mito_expanded_window <- plotting_window(c("V2", "V145"), nudge = 3)

mito_expanded_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "mito_expanded.pdf"), width = 3, height = 3)


# chrom expanded  -------------------------------------------------------

chrom_expanded_window <- plotting_window(c("V33", "V132"), nudge = 1)

chrom_expanded_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "chrom_expanded.pdf"), width = 3, height = 3)


# misc expanded  -------------------------------------------------------

misc_expanded_window <- plotting_window(c("V67", "V60"), nudge = 0.6)

misc_expanded_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "misc_expanded.pdf"), width = 3, height = 3)



# cyto expanded  -------------------------------------------------------

cyto_expanded_window <- plotting_window(c("V176", "V172"), nudge = 0.5)

cyto_expanded_window %>%
  plot_fn_panel() +
  theme(legend.position = "none") +
  ggsave(file.path(out_path, "cyto_expanded.pdf"), width = 3, height = 3)



# SHOC2 factors -----------------------------------------------------------
kras <- "V138"
nras <- "V209"
egfr <- "V162"
fgfr <- "V8"

#KRAS and NRAS
ras_window <- plotting_window(c(kras, nras), nudge = 0.1)

plot_loadings_panel(ras_window, kras) +
  ggsave(file.path(out_path, "shoc2_kras.pdf"), width =1.5, height = 1.5)
plot_loadings_panel(ras_window, nras)  +
  ggsave(file.path(out_path, "shoc2_nras.pdf"), width =1.5, height = 1.5)

#EGFR and FGFR
gfr_window  <- plotting_window(c(egfr, fgfr), nudge = 0.1)


plot_loadings_panel(gfr_window, egfr) +
  ggsave(file.path(out_path, "shoc2_egfr.pdf"), width =1.5, height = 1.5)
plot_loadings_panel(gfr_window, fgfr) +
  ggsave(file.path(out_path, "shoc2_fgfr.pdf"), width =1.5, height = 1.5)


param_in <- 2.2
#Whole embedding plots
plot_loadings_panel(loc_plot_df, kras, gene_labs = F) +
  ggsave(file.path(out_path, "global_kras.pdf"), width =param_in, height = param_in)

plot_loadings_panel(loc_plot_df, nras, gene_labs = F) +
  ggsave(file.path(out_path, "global_nras.pdf"), width =param_in, height = param_in)


plot_loadings_panel(loc_plot_df, egfr, gene_labs = F) +
  ggsave(file.path(out_path, "global_egfr.pdf"), width =param_in, height = param_in)

plot_loadings_panel(loc_plot_df, fgfr, gene_labs = F) +
  ggsave(file.path(out_path, "global_fgfr.pdf"), width =param_in, height = param_in)
