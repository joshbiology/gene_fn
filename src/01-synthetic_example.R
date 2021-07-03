#synthetic data for factorization

library(ProjectTemplate); load.project()
library(pheatmap)
library(sn)

out_path <- "./output/01-synthetic_example"

# Functions ---------------------------------------------------------------


generate_skewed_dictionary <- function(n = 25, rho = 0.2) {
  #https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  n    # length of vector
  rho   # desired correlation = cos(angle)
  theta <- acos(rho)             # corresponding angle
  x1    <- -dsn(seq(-0.5, 1, length.out = num_cell_lines) %>% sample, xi = 0.1, omega = 0.3, alpha = 5)       # fixed given data
  x2    <- -dsn(seq(-0.5, 1, length.out = num_cell_lines) %>% sample, xi = 0.1, omega = 0.3, alpha = 5)
  # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  cor(x1, x)                                    # check correlation = rho
  return(cbind(x1, x) %>% scale(center= T, scale = F))
}

# Set up synthetic cell lines and genes  ---------------------------------------------------
num_cell_lines <- 25
cl_names <- paste("Cell_Line", 1:num_cell_lines, sep = "_")

num_genes <- 100
gene_names <- paste("Gene", 1:num_genes, sep = "_")


gene_df <- tibble(Gene = gene_names,
                  Function = c(rep("1", 40), rep("2", 40), rep("Both", 20)) %>% 
                    factor(levels = c("1", "2", "Both")))


seq(-0.1, 1.1, length.out = num_cell_lines)



# Prepare data ------------------------------------------------------

#Option 1: For reproducing the figures:
dict_mat <- read_csv("./data/interim/2021-1-28_synthetic_dict.csv") %>% as.matrix()

#Option 2: For generating new data to try:
#dict_mat <- generate_skewed_dictionary(num_cell_lines, rho = 0.25) %>% normalize_cols()

#Plot fitness scatter
dict_mat %>% as_tibble %>% 
  mutate(Diff = (V1-V2)) %>%   
  ggplot(aes(V1, V2, color = Diff)) +
  geom_point(shape = 18, size = 2) +
  scale_color_gradient2(high = ("#EA89B8"), mid = "gray90", low = ("#4278BC"), limits = c(-0.2, 0.2), oob=scales::squish) +
  labs(color = "", 
       x = "Function 1 Dependency",
       y = "Function 2 Dependency")

# Generate synthetic data; plot heatmaps ----------------------------------

#Dictionary
pheatmap(dict_mat[(order(dict_mat[,1]-dict_mat[,2])),], cluster_cols = F, cluster_rows = F, colorRampPalette(c("red", "white"))(100),
         filename = "./output/figures/synthetic_dict.pdf", width = 2, height = 10)


code_mat <- rbind(c(rep(2, 40), rep(0, 40), rep(1, 20)),
                  c(rep(0, 40), rep(2, 40), rep(1, 20)))

reconstruct_mat <- dict_mat %*% (code_mat)

#Loadings
pheatmap(code_mat %>% set_colnames(gene_names), show_colnames = F, cluster_cols = F, cluster_rows = F,  colorRampPalette(c("white", "#00AEEF"))(100),
         annotation_col = gene_df %>% column_to_rownames("Gene"),
         annotation_colors = list(Function = c("1" = "#4278BC", "2" = "#EA89B8", "Both" = "#8575AD")), 
         filename =  file.path(out_path, "synthetic_g2f.pdf"), width = 10, height = 1, border_color = "gray80")

#Synthetic fitness data
pheatmap(dict_mat %*% (code_mat) %>% set_colnames(gene_names) %>% 
           magrittr::extract(order(dict_mat[,1]-dict_mat[,2]),), show_colnames = F, cluster_cols = F, cluster_rows = F, 
         annotation_col = gene_df %>% column_to_rownames("Gene"),
         annotation_colors = list(Function = c("1" = "#4278BC", "2" = "#EA89B8", "Both" = "#8575AD")),
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100),
         breaks = seq(-1, 1, length.out=101), 
         filename =  file.path(out_path, "synthetic_dep.pdf"), width = 10, height = 4)

#Noise
noise_mat <- rnorm(length(reconstruct_mat), sd = 0.3) %>% matrix(dim(reconstruct_mat))

final_mat <- (reconstruct_mat + noise_mat)

pheatmap(final_mat %>% set_colnames(gene_names) %>% 
           magrittr::extract(order(dict_mat[,1]-dict_mat[,2]),), show_colnames = F, cluster_cols = F, cluster_rows = F, 
         annotation_col = gene_df %>% column_to_rownames("Gene"),
         annotation_colors = list(Function = c("1" = "#4278BC", "2" = "#EA89B8", "Both" = "#8575AD")),
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100),
         breaks = seq(-1, 1, length.out=101), 
         filename =  file.path(out_path, "synthetic_dep+noise.pdf"), width = 10, height = 4)


# UMAP ---------------------------------------------------------
umap_out <- umap::umap(final_mat %>% t(), metric = "cosine", n_neighbors = 10, random_state = 123)



umap_out %>% 
  umap_to_df("Index") %>% 
  cbind(gene_df) %>% 
  ggplot(aes(V1, V2, color = Function)) +
  geom_point()


# Other factorizations ---------------------------------------------------------------------

pc <- prcomp(final_mat)

ggplot(dict_mat %>% as_tibble, aes(V1, V2)) +geom_point()

plot(pc$x[,1:2])

umap_out %>% 
  umap_to_df("Index") %>% 
  cbind(gene_df) %>% 
  cbind(pc$rotation[,1:2]) %>% 
  pivot_longer(names_to  = "Factor", values_to = "Loadings", starts_with("PC")) %>% 
  ggplot(aes(V1, V2, color = Loadings)) +
  geom_point() +
  facet_wrap(~Factor) +
  scale_color_viridis_c()

ic <- icasso_ica(final_mat, nbComp = 2)
plot(ic$S)

umap_out %>% 
  umap_to_df("Index") %>% 
  cbind(gene_df) %>% 
  cbind(ic$A %>% set_colnames(c("IC1", "IC2"))) %>% 
  pivot_longer(names_to  = "Factor", values_to = "Loadings", starts_with("IC")) %>% 
  ggplot(aes(V1, V2, color = Loadings)) +
  geom_point() +
  facet_wrap(~Factor) +
  scale_color_viridis_c()


k_m <- kmeans(final_mat %>% t(), centers = 2)

plot(k_m$centers %>% t())

umap_out %>% 
  umap_to_df("Index") %>% 
  cbind(gene_df) %>% 
  mutate(Cluster = factor(k_m$cluster)) %>% 
  ggplot(aes(V1, V2, color = Cluster)) +
  geom_point()



# Dictionary learning -----------------------------------------------------------------
#export_for_matlab(final_mat, "./output/synthetic_example.csv")
#Run DGRDL with matlab wrapper. I used the following code:
#/Applications/MATLAB_R2019b.app/bin/matlab -batch "DGRDL_wrapper('/Users/joshpan/gene_function_dev/output/synthetic_example.csv', 'K', 2, 'T', 2, 'filename_out', '/Users/joshpan/gene_function_dev/data/interim/matlab/synthetic/synthetic-K=2-T=2-n_neigh=1.mat', 'num_neighbor', 1);exit"

#Factorization output (patern accepts many files, but we only have one in this folder.)
mat_paths <- list.files("./data/interim/matlab/synthetic", pattern = "*.mat", full.names = T) %>% unlist()

grid_factorized <- map(mat_paths, function(x) {
  tmp <- R.matlab::readMat(x)
  tmp$X <- matrix(tmp$X, ncol = dim(tmp$X)[2])
  
  graphdl_to_factorized(import_graphdl(tmp), gene_names, cl_names)
})


params_df <- map_df(grid_factorized, function(x){
  tibble(K = attr(x, 'extras')$K %>% as.numeric,
         T_param = attr(x, 'extras')$T %>% as.numeric,
         Alpha = attr(x, 'extras')$alpha %>% as.numeric,
         Beta = attr(x, 'extras')$beta %>% as.numeric)
})

#Align to original function orderin (w/o loss of generality)
recovered_dictionary <- grid_factorized[[1]] %>% 
  get_cell_mat() %>% 
  magrittr::extract(,c(2,1))

recovered_loadings <-  grid_factorized[[1]] %>% 
  get_gene_mat() %>% 
  magrittr::extract(,c(2,1))


#Heatmaps

pheatmap(recovered_dictionary[(order(dict_mat[,1]-dict_mat[,2])),], cluster_cols = F, cluster_rows = F, show_rownames = F,colorRampPalette(c("red", "white"))(100),
         filename = file.path(out_path, "synthetic_graphdl_dict.pdf"), width = 2, height = 10, breaks = seq(-0.45, 0.15, length.out=101))

pheatmap(t(recovered_loadings), show_colnames = F, cluster_cols = F, cluster_rows = F,  colorRampPalette(c("white", "#00AEEF"))(100),
         annotation_col = gene_df %>% column_to_rownames("Gene"),
         annotation_colors = list(Function = c("1" = "#4278BC", "2" = "#EA89B8", "Both" = "#8575AD")), display_numbers = F,
         filename =  file.path(out_path, "synthetic_graphdl_g2f.pdf"), width = 25, height = 1, border_color = "gray80")


# Cell state scatters -----------------------------------------------------
tmp_names <- c("Ground_Truth", "PCA", "ICA", "K_Means", "Webster")

list(Ground_Truth = dict_mat,
     PCA = pc$x[,1:2] %*% diag(c(-1, -1)),
     ICA = ic$S,
     K_Means = k_m$centers %>% t() %>% magrittr::extract(,c(1,2)), #Depending on the stochasticity of k-means, you might have to flip columns
     Webster = recovered_dictionary) %>% 
  map(function(x) x %>% 
        scale() %>% 
        set_colnames(c("V1", "V2")) %>% 
        as_tibble() %>% 
        mutate(Name = cl_names)) %>% 
  enframe("Method") %>% 
  unnest(value) %>% 
  mutate(Method = factor(Method, levels = tmp_names)) %>% 
  left_join(tibble(Name = cl_names,
                   Diff = dict_mat[,1] - dict_mat[,2])) %>% 
  ggplot(aes(V1, V2, color = Diff)) +
  facet_wrap(~Method, nrow = 1) +
  geom_point(shape = 18, size = 4) +
  scale_color_gradient2(high = ("#EA89B8"), mid = "gray90", low = ("#4278BC"), limits = c(-0.2, 0.2), oob=scales::squish) +
  labs(color = "", 
       x = "Function 1 Dependency",
       y = "Function 2 Dependency") +
  theme_minimal() +
  ggsave(file.path(out_path, "recovering_cell_state_synthetic.pdf"), width = 15, height = 3.5)


# UMAP scatters -----------------------------------------------------------

list(Ground_Truth = code_mat %>% t(),
     K_Means = (model.matrix(~0+k_m$cluster %>% as.factor()) *2 ) %>% matrix(ncol = 2) %>% magrittr::extract(,c(2,1)),
     Webster = recovered_loadings) %>% 
  map(~as_tibble(.) %>% mutate(Gene = gene_names)) %>% 
  enframe() %>% 
  unnest(value) %>% 
  set_colnames(c("Dataset", "Function_1", "Function_2", "Gene")) %>% 
  left_join(
    umap_out %>% 
      umap_to_df("Index") %>% 
      cbind(gene_df)) %>% 
  pivot_longer(names_to  = "Factor", values_to = "Loadings", c("Function_1", "Function_2")) %>% 
  ggplot(aes(V2, V1, color = Loadings)) +
  geom_point(size = 2.5) +
  facet_grid(Factor~Dataset) +
  scale_color_gradient2(low = "gray90", mid = ("#a184d1"), high = "#00AEEF", limits = c(0, 2), midpoint = 1, oob = scales::squish)  +
  theme_void() +
  ggsave(file.path(out_path, "recovering_gene_manifold_synthetic.pdf"), width = 15, height = 3.5, device= cairo_pdf)
