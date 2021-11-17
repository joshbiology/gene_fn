library(ProjectTemplate); load.project()
require("caret")

out_path <- file.path(".", "output", "09-projection")
create_output_folder(out_path)

# Functions ---------------------------------------------------------------

impute_prism_na <- function(pri) {
  require(FastImputation)
  X =  pri  %>% t() %>% as.data.frame(); X$id = rownames(X)
  patterns <- TrainFastImputation(X,idvars="id")
  X <- FastImputation(X, patterns, verbose = T)
  rownames(X) = X$id; X$id = NULL
  X = as.matrix(X) %>% t()
  rownames(X) <- rownames(pri)
  colnames(X) <- colnames(pri)
  return(X)
}


#Modes package deprecated from CRAN
#Source code to allow modes::bimodality_coefficient to run
#https://github.com/sathish-deevi/modes-Package/blob/d4a95a49dc86d434fdfd9b9d541d78c58dd4a45a/src/R/Utility_functions.R
skewness<-function(x, finite=TRUE){
  n=length(x)
  S=(1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  if(finite==FALSE){
    S=S
  }else{
    S=S*(sqrt(n*(n-1)))/(n-2)
  }
  return(S)	
}


kurtosis<-function(x, finite){
  n=length(x)
  K=(1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2) - 3
  if(finite==FALSE){
    K=K
  }
  else{
    K=((n-1)*((n+1)*K - 3*(n-1))/((n-2)*(n-3))) +3
  }
  return(K)	
}



#https://github.com/sathish-deevi/modes-Package/blob/d4a95a49dc86d434fdfd9b9d541d78c58dd4a45a/src/R/Nonparametric_functions.R
bimodality_coefficient<-function(x, finite=TRUE,...){
  if(finite==TRUE){
    G=skewness(x,finite)
    sample.excess.kurtosis=kurtosis(x,finite)
    K=sample.excess.kurtosis
    n=length(x)
    B=((G^2)+1)/(K+ ((3*((n-1)^2))/((n-2)*(n-3))))
  }
  else{
    G=skewness(x,FALSE)
    K=kurtosis(x,FALSE)
    B=((G^2)+1)/(K)
  }
  return(B)
}

# Data -------------------------------------------------------------
load("./cache/avana_19q4_webster.RData")
source("./munge/webster_depmap.R")
source("./munge/prism.R")


# Preprocess --------------------------------------------------------------


#Remove cell line center (trimmed means)
prism_means <- prism_lfc %>% aaply(1, function(x) mean(x, trim = 0.1, na.rm = T)) %>% set_names(rownames(prism_lfc))
prism_lfc_centered <- prism_lfc - prism_means

#Apply trimmed means from primary data to secondary data
prism_secondary_means <- prism_means[rownames(prism_secondary_lfc)]
prism_secondary_centered <- prism_secondary_lfc - prism_secondary_means


# Primary data exploration ------------------------------------------------

#Curate primary screen compounds
prism_primary_compounds <- prism_umap_annot %>%
  dplyr::filter(group != "OTHER") %>%
  pull(column_name) %>% unique()

pri_tmp <- prism_lfc_centered[,prism_primary_compounds]

pri.bc = apply(pri_tmp,  2, function(x) bimodality_coefficient(x[is.finite(x)]))
pri.ns = apply(pri_tmp, 2, function(x) sum(x  < log2(.3), na.rm  = T))
pri.cl = apply(pri_tmp, 1, function(x) sum(is.na(x), na.rm  = T))
pri_tmp = pri_tmp[pri.cl < nrow(pri_tmp)*0.2, (pri.bc > 0.2) &  (pri.ns > 4)]

#Impute
prism_primary_imputed <- impute_prism_na(pri_tmp)

ProjectTemplate::cache("prism_primary_imputed")

#For reproducibility:
#load("./cache/prism_primary_imputed.Rdata")

#Calculate statitics
primary_compound_stats_df <- tibble(column_name = colnames(prism_primary_imputed),
                                    Var = matrixStats::colVars(prism_primary_imputed),
                                    Bimodality = map_dbl(column_name, function(x) bimodality_coefficient(prism_primary_imputed[,x])),
                                    Median = matrixStats::colMedians(prism_primary_imputed))



# Plot compound landscape -------------------------------------------------

umap_compounds_orig <- umap::umap(prism_primary_imputed %>% t(), metric = "cosine", n_neighbors = 5, random_state = 1)

g1 <- umap_compounds_orig$layout %>%
  as_tibble() %>%
  dplyr::mutate(column_name = colnames(prism_primary_imputed)) %>%
  left_join(primary_compound_stats_df) %>%
  left_join(prism_umap_annot %>%
              dplyr::filter(group != "OTHER")) %>%
  ggplot(aes(V1, V2, color = group, text = paste(group, name, Median))) +
  geom_point()

#plotly::ggplotly(g1)

umap_compounds_orig$layout %>%
  as_tibble() %>%
  dplyr::mutate(column_name = colnames(prism_primary_imputed)) %>%
  left_join(prism_umap_annot %>%
              dplyr::filter(group != "OTHER")) %>%
  ggplot(aes(V1, V2, color = group)) +
  geom_point() +
  facet_wrap(~group) +
  theme_minimal() +
  theme(legend.position = "NA")

#Plot median cell fitness
umap_compounds_orig$layout %>%
  as_tibble() %>%
  dplyr::mutate(column_name = colnames(prism_primary_imputed)) %>%
  left_join(primary_compound_stats_df) %>%
  ggplot(aes(V1, V2, color = Median)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

umap_compounds_orig$layout %>%
  as_tibble() %>%
  dplyr::mutate(column_name = colnames(prism_primary_imputed)) %>%
  left_join(primary_compound_stats_df) %>%
  ggplot(aes(V1, V2, color = Median < -1.5)) +
  geom_point()

primary_compound_filtered_meta <- primary_compound_stats_df %>%
  dplyr::filter(Median > -1.5) %>%
  left_join(prism_umap_annot %>%
              dplyr::filter(group != "OTHER")) %>%
  add_count(group, name = "Num_Compounds") %>%
  dplyr::filter(Num_Compounds >= 5) %>%
  arrange(group, name)

prism_primary_compounds_omp <- unique(primary_compound_filtered_meta$column_name)


ProjectTemplate::cache("primary_compound_filtered_meta")

# Prepare primary data for OMP ------------------------
#Orthogonal matching pursuit requires 1. a pre-learned dictionary. 2. a set of new "signals" to model in terms of dict elements
primary_common_cl <- get_cell_mat(webster_depmap) %>% rownames() %>% intersect(rownames(prism_primary_imputed))

prism_primary_signal_input <- prism_primary_imputed[primary_common_cl,prism_primary_compounds_omp]

prism_primary_dict_input <- get_cell_mat(webster_depmap)[primary_common_cl,]


# Quick qc ----------------------------------------------------------------

umap_filtered <- umap::umap(prism_primary_signal_input %>% t(), metric = "cosine", n_neighbors = 5, random_state = 1)

g1 <- umap_filtered$layout %>%
  as_tibble() %>%
  dplyr::mutate(column_name = prism_primary_compounds_omp) %>%
  left_join(primary_compound_stats_df) %>%
  left_join(prism_umap_annot %>%
              dplyr::filter(group != "OTHER")) %>%
  ggplot(aes(V1, V2, color = group, text = paste(group, name, Median))) +
  geom_point()

#plotly::ggplotly(g1)

# Prepare orthogonal matching pursuit -------------------------------------
#This code block prepares a bash script to call MATLAB and perform orthogonal matching pursuit.
#For the repo, I'm including a flag to skip this step and read from the stored output directly.

run_local_omp <- F

if (run_local_omp) {
  
  prism_primary_dict_path <- file.path(out_path, "prism_primary_dict_input.csv")
  prism_primary_signal_path <- file.path(out_path, "prism_primary_signal_input.csv")
  prism_primary_omp_file <- file.path(out_path, "prism_primary_omp_output.csv")
  
  #Export raw data
  export_for_matlab(prism_primary_dict_input, out_file =prism_primary_dict_path)
  export_for_matlab(prism_primary_signal_input, out_file = prism_primary_signal_path)
  
  #T_param
  T_param <- 4
  
  #Modify to your local matlab path
  sys_command <- sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
                         prism_primary_dict_path,
                         prism_primary_signal_path,
                         prism_primary_omp_file,
                         T_param)
  
  cat(sys_command, sep = '\n',file=file.path(".", "cluster_scripts", paste("prism_primary_omp","txt", sep = ".")))
  
  #Run the script. You can also do this in the terminal
  system("bash ./cluster_scripts/prism_primary_omp.txt")
  
}



# Process PRISM primary OMP -----------------------------------------------
#As mentioned, the repo will contain the OMP results in the interim folder.
#This file is also available on the FigShare:
#https://figshare.com/articles/dataset/Webster_Supplemental_Output/14963561?file=29073813

prism_primary_omp <- read_tsv(file.path(".", "data", "interim", "prism_primary_omp.tsv"), col_names = T) %>%
  column_to_rownames("column_name") %>%
  as.matrix()

ProjectTemplate::cache("prism_primary_omp")

recon_prism <- tcrossprod(prism_primary_dict_input, prism_primary_omp)

proj_results <- tibble(column_name = colnames(recon_prism),
                       Pearson = map_dbl(column_name, function(x) cor(recon_prism[,x], prism_primary_signal_input[,x])))

ProjectTemplate::cache("proj_results")

proj_results %>%
  left_join(primary_compound_filtered_meta) %>%
  ggplot(aes(reorder(group, Pearson, FUN = median), Pearson, fill = reorder(group, Pearson, FUN = median))) +
  geom_boxplot() + coord_flip() +
  theme(legend.position = "none")

ggsave(file.path(out_path, "primary_pearson_recon.pdf"), width = 7, height = 4, device = cairo_pdf)

proj_results %>%
  left_join(primary_compound_filtered_meta) %>%
  ggplot(aes(Median, Pearson, color = reorder(group, Pearson, FUN = median))) +
  geom_point()

ggsave(file.path(out_path, "primary_pearson_recon_scatter.pdf"), width = 8, height = 4, device = cairo_pdf)

# Plot loadings -----------------------------------------------------------
plot_moa_loadings <- function(selected_group) {
  selected_compounds <- primary_compound_filtered_meta %>% dplyr::filter(group == selected_group) %>% pull(column_name)

  tmp_omp <- prism_primary_omp[selected_compounds,]/sqrt(365)


  selected_fns <- colSums(tmp_omp != 0) %>%
    enframe() %>%
    arrange(-value) %>%
    dplyr::filter(value > 2) %>%
    slice(1:10) %>%
    pull(name)

  show_numbers <- length(selected_compounds) < 40

  annot <- proj_results %>% left_join(primary_compound_filtered_meta) %>%
    group_by(column_name) %>%
    slice(1) %>%
    ungroup() %>% column_to_rownames("column_name") %>% dplyr::select(dose, name, screen_id, Pearson, Median)

  pheatmap::pheatmap(tmp_omp[,selected_fns], show_rownames = F, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", fontsize_col = 7.5,
                     annotation_row = annot,
                     color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), main = selected_group, display_numbers = show_numbers, number_format = "%.2f",
                     breaks = seq(-0.5, 0.5, length.out=101), filename = file.path(out_path, paste(selected_group, "loading_heatmap.pdf", sep = "_")))

}

walk(primary_compound_filtered_meta$group %>% unique() %>% setdiff(c("HIF PROLYL HYDROXYLASE", "RET TYROSINE KINASE INHIBITOR")), plot_moa_loadings)


# Project imputed profiles ------------------------------------------------
#select successful moa's based on pearson
successful_moas <- proj_results %>%
  left_join(primary_compound_filtered_meta) %>%
  group_by(moa_umap) %>%
  summarize(Median_Pearson = median(Pearson)) %>%
  dplyr::mutate(Successful = Median_Pearson > 0.33)

#select successful compounds's based on pearson
successful_cmpds <-proj_results %>%
  left_join(primary_compound_filtered_meta) %>%
  dplyr::filter(moa_umap %in% (successful_moas %>% dplyr::filter(Successful) %>% pull(moa_umap))) %>%
  dplyr::filter(Pearson >= .33) %>%
  pull(column_name) %>%
  unique()

#Impute onto full dictionary
compound_imputed <- tcrossprod(get_cell_mat(webster_depmap), prism_primary_omp)

#Project onto original embedding
load("./cache/depmap_umap.RData")
compound_layout <- predict(depmap_umap, compound_imputed[, successful_cmpds] %>% t())

#Visualize
impute_df <- depmap_umap$layout %>% rbind(compound_layout) %>%
  as_tibble(rownames = "Name") %>%
  dplyr::mutate(Type = c(rep("Function", webster_depmap$rank), rep("Gene", dim(avana_19q4_webster)[2]), rep("Compound", length(successful_cmpds) ))) %>%
  dplyr::mutate(Gene_Name = convert_genes(Name, from = "entrez_id", "symbol")) %>%
  rename(X = V1, Y = V2) %>%
  dplyr::mutate(Y = -Y)  %>%
  left_join(primary_compound_filtered_meta, by = c("Name" = "column_name")) %>%
  dplyr::mutate(moa_umap  =factor(moa_umap))

ProjectTemplate::cache("impute_df")

g_combo <- impute_df %>%
  ggplot(aes(X, Y, label = paste(Name, name))) +
  geom_point(data = subset(impute_df, Type == "Gene"), size = 0.5, alpha = 0.2, color = "gray70") +
  geom_point(data = subset(impute_df, Type == "Function"),size = 1, shape = 17, alpha = 0.5) +
  geom_point(data = subset(impute_df, Type == "Compound"),size = 2, shape = 18,aes(color = moa_umap)) +
  theme_void()

#plotly::ggplotly(g_combo)

g_combo

ggsave(file.path(out_path, "compound_embedding.pdf"), width = 5.5, height = 4, device = cairo_pdf)




# Focus insets ------------------------------------------------------------

#adapted from 08-landscape.R
plotting_window <- function(fn_names, nudge = 0.3) {

  tmp <- impute_df %>% dplyr::filter(Name %in% fn_names)
  window_x <- c(min(tmp$X)-nudge, max(tmp$X)+nudge)
  window_y <- c(min(tmp$Y)-nudge, max(tmp$Y)+nudge)
  impute_df %>%
    dplyr::filter(X > window_x[1], X < window_x[2], Y > window_y[1], Y < window_y[2])
}

plot_fn_panel <- function(df) {
  df %>%
    ggplot(aes(X, Y, label = Name)) +
    geom_point(data = subset(df, Type == "Gene"), size = 0.5, alpha = 0.5, color = "gray70") +
    geom_point(data = subset(df, Type == "Function"),size = 1, shape = 17, alpha = 0.5) +
    geom_point(data = subset(df, Type == "Compound"),size = 2, shape = 18,aes(color = moa_umap)) +
    geom_text(data = df %>% dplyr::filter(Type == "Function"),  nudge_y = .05) +
    scale_colour_hue(drop=FALSE) +
    scale_alpha(limits = c(0, 0.6), range = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none")
}


plotting_window(c("V15"), 0.1) %>% plot_fn_panel

ggsave(file.path(out_path, "braf.pdf"), width = 1.5, height = 1.5, device = cairo_pdf)

plotting_window(c("V151", "BRD-K69726342-238-02-4::2.5::HTS", "BRD-K94441233-001-13-0::2.6::HTS", "BRD-K94441233-001-17-1::2.5::HTS", "BRD-K22134346-001-24-9::2.42::HTS"), 0.1) %>% plot_fn_panel

ggsave(file.path(out_path, "statin.pdf"), width = 1.5, height = 1.5, device = cairo_pdf)

plotting_window(c("V109", "BRD-K44665581-001-01-8::2.5::HTS"), 0.3) %>% plot_fn_panel

ggsave(file.path(out_path, "brom.pdf"), width = 1.5, height = 1.5, device = cairo_pdf)

plotting_window(c("V82"), 0.15) %>% plot_fn_panel

ggsave(file.path(out_path, "egfr.pdf"), width = 1.5, height = 1.5, device = cairo_pdf)


# Compound to function heatmaps -------------------------------------------

plot_compound_loadings <- function(index, max_loading = F) {
  cmpd_df <- prism_primary_omp[,index] %>%
    enframe("column_name", "Loadings")

  cmpd_df <- cmpd_df %>%
    dplyr::mutate(Rank = order(order(Loadings, decreasing = T))) %>%
    arrange(desc(Loadings)) %>%
    dplyr::filter(Rank <= 10) %>%
    left_join(prism_umap_annot %>% group_by(column_name) %>% dplyr::filter(row_number()==1) %>% ungroup)

  loading_mat <- prism_primary_omp[cmpd_df$column_name,index] %>% matrix(ncol = 1) %>%
    set_rownames(cmpd_df$column_name)

  loading_mat <- loading_mat/sqrt(367) #converting euclidean to stdev

  if (max_loading == F) {
    max_loading <- max(abs(loading_mat))
  }

  g2 <- pheatmap::pheatmap(loading_mat, cluster_cols = F, cluster_rows = F, show_colnames = F,
                           annotation_row = cmpd_df %>% dplyr::select(column_name, moa_umap) %>% column_to_rownames("column_name"),
                           labels_row = cmpd_df$name,
                           filename =  file.path(out_path, paste("cmpd_loading_", index, ".pdf", sep = "")),
                           cellwidth =20, cellheight = 20, color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100),
                           breaks = seq(-max_loading, max_loading, length.out=101))
}

plot_compound_loadings(151)
plot_compound_loadings(15)
plot_compound_loadings(82)
plot_compound_loadings(109)

# Secondary screen --------------------------------------------------------

secondary_meta <- prism_secondary_meta %>% dplyr::filter(moa %in% (successful_moas %>% dplyr::filter(Successful) %>% pull(moa_umap))) %>% #filter(!(grepl("PR300", compound_plate))) %>%
  arrange(moa, name, screen_id, dose)

ProjectTemplate::cache('secondary_meta')

prism_secondary_tmp <- prism_secondary_lfc[,secondary_meta$column_name]


# Preprocess prism secondary screening data -------------------------------
#Impute remaining data. Takes a few mins
prism_secondary_imputed <- impute_prism_na(prism_secondary_tmp)
ProjectTemplate::cache("prism_secondary_imputed")

#For reproducibility:
#load('./cache/prism_secondary_imputed.RData')

# Prepare PRISM matrix with CRISPR dictionary ------------------------------------
secondary_common_cl <- intersect(primary_common_cl, rownames(prism_secondary_imputed))

prism_omp_input <- prism_secondary_imputed[secondary_common_cl,]

dict_input <- get_cell_mat(webster_depmap)[secondary_common_cl,]

# Prepare orthogonal matching pursuit -------------------------------------
#As with the primary OMP, this code block connects to MATLAB. For the repo,
#I've shortcuted this step and provided the output below.

run_secondary_omp <- F

if (run_secondary_omp) {
  
  dict_path <- file.path(out_path, "prism_secondary_dict_input.csv")
  signal_path <- file.path(out_path, "prism_secondary_signal_input.csv")
  omp_file <- file.path(out_path, "prism_secondary_omp_output.csv")
  
  #Export raw data
  export_for_matlab(dict_input, out_file =dict_path)
  export_for_matlab(prism_omp_input, out_file = signal_path)
  
  #T_param
  T_param <- 4
  
  #Modify to your local matlab path
  sys_command <- sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
                         dict_path,
                         signal_path,
                         omp_file,
                         T_param)
  
  cat(sys_command, sep = '\n',file=file.path(".", "cluster_scripts", paste("prism_secondary_omp","txt", sep = ".")))
  
  #Run the script. You can also do this in the terminal
  system("bash ./cluster_scripts/prism_secondary_omp.txt")
  
}


# Results -----------------------------------------------------------------
#The result of the above code block is stored in gene_fn/data/interim, 
#and can also be downloaded here:
#https://figshare.com/articles/dataset/Webster_Supplemental_Output/14963561?file=29073804

prism_secondary_omp <- read_tsv(file.path(".", "data", "interim", "prism_secondary_omp.tsv"), col_names = T)

ProjectTemplate::cache("prism_secondary_omp")

recon_prism <- tcrossprod(dict_input, prism_secondary_omp %>% select(-column_name) %>% as.matrix())

proj_results_secondary <- secondary_meta %>% 
  select(column_name, compound_plate) %>% 
  mutate(Pearson = map_dbl(1:ncol(recon_prism), function(x) cor(recon_prism[,x], prism_omp_input[,x])),
         Median = matrixStats::colMedians(prism_omp_input))

ProjectTemplate::cache("proj_results_secondary")

# Plot secondary loadings -----------------------------------------------------------
plot_secondary_moa_loadings <- function(selected_group) {
  #In secondary_meta, the column_name and compound_plate form a unique ID for each compound. 
  
  tmp_meta <- secondary_meta %>% 
    unite("ID", c("column_name", "compound_plate"))
  
  selected_compounds <- tmp_meta %>% 
    dplyr::filter(moa == selected_group) %>% 
    pull(ID)

  tmp_omp <- prism_secondary_omp[tmp_meta$ID %in% selected_compounds,] %>% 
    select(-column_name) %>% 
    as.matrix() %>% 
    set_rownames(selected_compounds)


  selected_fns <- colSums(tmp_omp != 0) %>%
    enframe() %>%
    arrange(-value) %>%
    dplyr::filter(value > 2) %>%
    slice(1:10) %>%
    pull(name)

  show_numbers <- length(selected_compounds) < 40

  annot <- proj_results_secondary %>% 
    unite("ID", c("column_name", "compound_plate")) %>%
    left_join(tmp_meta) %>%
    group_by(ID) %>%
    slice(1) %>%
    ungroup() %>% column_to_rownames("ID") %>% dplyr::select(dose, name, screen_id, Pearson, Median)

  pheatmap::pheatmap(tmp_omp[,selected_fns], show_rownames = F, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", fontsize_col = 7.5,
                     annotation_row = annot,
                     color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), main = selected_group, display_numbers = show_numbers, number_format = "%.0f",
                     breaks = seq(-10, 10, length.out=101), filename = file.path(out_path, paste("secondary", selected_group,  "loading_heatmap.pdf", sep = "_")))

}

walk(successful_moas %>% dplyr::filter(Successful) %>% pull(moa_umap), plot_secondary_moa_loadings)


# Plot individual loadings per dose ---------------------------------------

plot_loadings_by_dose <- function(moa_focus, fn, number_of_facets = 3) {
  #sqrt(number of data points) converts OMP coefficients to "loadings"
  
  cmpd_loadings <- prism_secondary_omp %>% select(-column_name) %>% as.matrix() %>% "/"(325)
  tmp <- proj_results_secondary %>%
    cbind(cmpd_loadings %>% as_tibble()) %>%
    left_join(secondary_meta)

  high_loaded_names <- tmp %>%
    dplyr::filter(moa == moa_focus) %>%
    group_by(name) %>%
    rename(Focus = sym(fn)) %>%
    summarize(Max_Loading = max(Focus)) %>%
    arrange(desc(Max_Loading)) %>%
    slice(1:number_of_facets) %>%
    pull(name)

  dose_df <- tmp %>%
    dplyr::filter(name %in% high_loaded_names) %>%
    pivot_longer(names_to = "Mode", values_to = "Value", c(sym(fn), "Median")) %>%
    dplyr::mutate(Mode = factor(Mode, levels = c(sym(fn), "Median")))



  dose_df %>%
    ggplot(aes(log(dose,2), Value, color)) +
    geom_point(aes(color = Mode), size = 0.5) +
    geom_col(width = 0.3, aes(fill = Mode)) +
    facet_grid(Mode~name, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Log2 Dose (uM)") +
    scale_color_manual(values=c("#2aace2", "red")) +
    scale_fill_manual(values=c("#2aace2", "red")) +
    theme(legend.position = "NA")
}


plot_loadings_by_dose("AKT inhibitor", "V61")

plot_loadings_by_dose("RAF inhibitor", "V15")

plot_loadings_by_dose("HMGCR inhibitor", "V151", 4)

ggsave(file.path(out_path, "hmgcr_secondary.pdf"), width = 7, height = 3, device = cairo_pdf)

plot_loadings_by_dose("MDM inhibitor", "V5", 2)

ggsave(file.path(out_path, "mdm_secondary.pdf"), width = 3.5, height = 3, device = cairo_pdf)

plot_loadings_by_dose("MEK inhibitor", "V138")

plot_loadings_by_dose("MEK inhibitor", "V15")

plot_loadings_by_dose("MEK inhibitor", "V82", 4)

ggsave(file.path(out_path, "meki_secondary.pdf"), width = 7, height = 3, device = cairo_pdf)



plot_loadings_by_dose("EGFR inhibitor", "V82", 4)

ggsave(file.path(out_path, "egfr_secondary.pdf"), width = 7, height = 3, device = cairo_pdf)

plot_loadings_by_dose("bromodomain inhibitor", "V109", 2)

ggsave(file.path(out_path, "bromo_secondary.pdf"), width = 3.5, height = 3, device = cairo_pdf)

