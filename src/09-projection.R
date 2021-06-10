library(ProjectTemplate); load.project()
require("caret")

out_path <- file.path(".", "output", "09-projection")

# Functions ---------------------------------------------------------------

impute_prism_na <- function(pri) {
  require(FastImputation)
  X =  pri  %>% t() %>% as.data.frame(); X$id = rownames(X)
  patterns <- TrainFastImputation(X,idvars="id")
  X <- FastImputation(X, patterns, verbose = T)
  rownames(X) = X$id; X$id = NULL
  X = as.matrix(X) %>% t()
}


# Data -------------------------------------------------------------
load("./cache/avana_19q4_webster.RData")
source("./munge/webster_depmap.R")
source("./munge/prism.R")


# Primary data exploration ------------------------------------------------


#Curate primary screen compounds
prism_primary_compounds <- prism_umap_annot %>% 
  filter(group != "OTHER") %>% 
  pull(column_name) %>% unique()

pri_tmp <- prism_lfc[,prism_primary_compounds]

#pri.bc = apply(pri_tmp,  2, function(x) modes::bimodality_coefficient(x[is.finite(x)]))
pri.ns = apply(pri_tmp, 2, function(x) sum(x  < log2(.3), na.rm  = T))
pri.cl = apply(pri_tmp, 1, function(x) sum(is.na(x), na.rm  = T))
pri_tmp = pri_tmp[pri.cl < nrow(pri_tmp)*0.2,   (pri.ns > 4)]#(pri.bc > 0.2) &

#Impute
prism_primary_imputed <- impute_prism_na(prism_lfc[,prism_primary_compounds])
colnames(prism_primary_imputed) <- prism_primary_compounds

#Calculate statitics
primary_compound_stats_df <- tibble(column_name = prism_primary_compounds,
       Var = matrixStats::colVars(prism_primary_imputed),
       Bimodality = map_dbl(column_name, function(x) modes::bimodality_coefficient(prism_primary_imputed[,x])),
       Median = matrixStats::colMedians(prism_primary_imputed))



# Plot compound landscape -------------------------------------------------

umap_compounds_orig <- umap::umap(prism_primary_imputed %>% t() %>% scale, metric = "cosine", n_neighbors = 5, random_state = 1)

g1 <- umap_compounds_orig$layout %>% 
  as_tibble() %>% 
  mutate(column_name = prism_primary_compounds) %>% 
  left_join(primary_compound_stats_df) %>% 
  left_join(prism_umap_annot %>% 
              filter(group != "OTHER")) %>% 
  ggplot(aes(V1, V2, color = group, text = paste(group, name, Median))) +
  geom_point()

plotly::ggplotly(g1)

umap_compounds_orig$layout %>% 
  as_tibble() %>% 
  mutate(column_name = prism_primary_compounds) %>% 
  left_join(prism_umap_annot %>% 
              filter(group != "OTHER")) %>% 
  ggplot(aes(V1, V2, color = group)) +
  geom_point() +
  facet_wrap(~group) +
  theme_minimal() +
  theme(legend.position = "NA") 

#Plot median cell fitness
umap_compounds_orig$layout %>% 
  as_tibble() %>% 
  mutate(column_name = prism_primary_compounds) %>% 
  left_join(primary_compound_stats_df) %>% 
  ggplot(aes(V1, V2, color = Median)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

umap_compounds_orig$layout %>% 
  as_tibble() %>% 
  mutate(column_name = prism_primary_compounds) %>% 
  left_join(primary_compound_stats_df) %>% 
  ggplot(aes(V1, V2, color = Median < -1.5)) +
  geom_point()

primary_compound_filtered_meta <- primary_compound_stats_df %>% 
  filter(Median > -1.5) %>% 
  left_join(prism_umap_annot %>% 
              filter(group != "OTHER")) %>% 
  add_count(group, name = "Num_Compounds") %>%  
  filter(Num_Compounds >= 5) %>% 
  arrange(group, name)

prism_primary_compounds_omp <- unique(primary_compound_filtered_meta$column_name)




# Prepare primary data for OMP ------------------------
#Orthogonal matching pursuit requires 1. a pre-learned dictionary. 2. a set of new "signals" to model in terms of dict elements
primary_common_cl <- get_cell_mat(webster_depmap) %>% rownames() %>% intersect(rownames(prism_primary_imputed))

prism_primary_signal_input <- prism_primary_imputed[primary_common_cl,prism_primary_compounds_omp] 

#scale cell line data 
#prism_primary_signal_input <- prism_primary_signal_input %>% t() %>% scale() %>% t()

prism_primary_dict_input <- get_cell_mat(webster_depmap)[primary_common_cl,]


# Quick qc ----------------------------------------------------------------

umap_filtered <- umap::umap(prism_primary_signal_input %>% t(), metric = "cosine", n_neighbors = 5, random_state = 1)

g1 <- umap_filtered$layout %>% 
  as_tibble() %>% 
  mutate(column_name = prism_primary_compounds_omp) %>% 
  left_join(primary_compound_stats_df) %>% 
  left_join(prism_umap_annot %>% 
              filter(group != "OTHER")) %>% 
  ggplot(aes(V1, V2, color = group, text = paste(group, name, Median))) +
  geom_point()

plotly::ggplotly(g1)

# Prepare orthogonal matching pursuit -------------------------------------

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



# Process PRISM primary OMP -----------------------------------------------

prism_primary_omp <- read_csv(prism_primary_omp_file, col_names = F) %>% 
  as.matrix() %>% 
  t() %>% 
  set_rownames(colnames(prism_primary_signal_input)) %>% 
  set_colnames(colnames(prism_primary_dict_input))

recon_prism <- tcrossprod(prism_primary_dict_input, prism_primary_omp)

proj_results <- tibble(column_name = colnames(recon_prism),
                       Pearson = map_dbl(column_name, function(x) cor(recon_prism[,x], prism_primary_signal_input[,x])))

proj_results %>% 
  left_join(primary_compound_filtered_meta) %>% 
  ggplot(aes(reorder(group, Pearson, FUN = median), Pearson)) +
  geom_boxplot() + coord_flip() +
  ggsave(file.path(out_path, "primary_pearson_recon.pdf"), width = 7, height = 7)


# Plot loadings -----------------------------------------------------------


plot_moa_loadings <- function(selected_group) {
  selected_compounds <- primary_compound_filtered_meta %>% dplyr::filter(group == selected_group) %>% pull(column_name)
  
  tmp_omp <- prism_primary_omp[selected_compounds,]
  
  
  selected_fns <- colSums(tmp_omp != 0) %>% 
    enframe() %>% 
    arrange(-value) %>% 
    filter(value > 2) %>% 
    slice(1:10) %>% 
    pull(name)
  
  show_numbers <- length(selected_compounds) < 40
  
  annot <- proj_results %>% left_join(primary_compound_filtered_meta) %>% 
    group_by(column_name) %>% 
    slice(1) %>% 
    ungroup() %>% column_to_rownames("column_name") %>% select(dose, name, screen_id, Pearson, Median)
  
  pheatmap::pheatmap(tmp_omp[,selected_fns], show_rownames = F, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", fontsize_col = 7.5,
                     annotation_row = annot,
                     color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), main = selected_group, display_numbers = show_numbers, number_format = "%.0f",
                     breaks = seq(-10, 10, length.out=101), filename = file.path(out_path, paste(selected_group, "loading_heatmap.pdf", sep = "_")))
  
}

walk(primary_compound_filtered_meta$group %>% unique() %>% setdiff(c("HIF PROLYL HYDROXYLASE", "RET TYROSINE KINASE INHIBITOR")), plot_moa_loadings)


plot_moa_loadings("RAF INHIBITOR")


# Project imputed profiles ------------------------------------------------
successful_moas <- proj_results %>% 
  left_join(primary_compound_filtered_meta) %>% 
  group_by(moa_umap) %>% 
  summarize(Median_Pearson = median(Pearson)) %>% 
  mutate(Successful = Median_Pearson > 0.33)

successful_cmpds <-proj_results %>% 
  left_join(primary_compound_filtered_meta) %>% 
  filter(moa_umap %in% (successful_moas %>% filter(Successful) %>% pull(moa_umap))) %>%
  pull(column_name) %>% 
  unique()

compound_imputed <- tcrossprod(get_cell_mat(webster_depmap), prism_primary_omp)

load("./cache/depmap_umap.RData")

compound_layout <- predict(depmap_umap, compound_imputed[, successful_cmpds] %>% t())

impute_df <- depmap_umap$layout %>% rbind(compound_layout) %>% 
  as_tibble(rownames = "Name") %>% 
  mutate(Type = c(rep("Factor", webster_depmap$rank), rep("Gene", dim(avana_19q4_webster)[2]), rep("Compound", length(successful_cmpds) ))) %>% 
  mutate(Gene_Name = convert_genes(Name, from = "entrez_id", "symbol")) %>% 
  rename(X = V1, Y = V2)


g_combo <- impute_df %>% 
  ggplot(aes(X, Y)) + 
  geom_point(data = subset(impute_df, Type == "Gene"), size = 0.5, alpha = 0.2, color = "gray70") +
  geom_point(data = subset(impute_df, Type == "Factor"),size = 1, shape = 17, alpha = 0.5) +
  geom_point(data = subset(impute_df, Type == "Compound") %>% left_join(primary_compound_filtered_meta, by = c("Name" = "column_name")),size = 2, shape = 18,aes(color = moa_umap)) +
  theme_void()

plotly::ggplotly(g_combo)

g_combo +
  ggsave("./output/primary_compound_embedding_combined.pdf", width = 9, height = 7, device = cairo_pdf)


plotly::ggplotly(g_combo)




# Secondary screen --------------------------------------------------------



# Curate compound MOA's ---------------------------------------------------

secondary_meta <- prism_secondary_meta %>% filter(moa %in%  (successful_moas %>% filter(Successful) %>% pull(moa_umap))) %>% filter(!(grepl("PR300", compound_plate))) %>% 
  arrange(moa, name, screen_id, dose)

prism_tmp <- prism_secondary_lfc[,secondary_meta$column_name]


# Preprocess prism secondary screening data -------------------------------
#Cell lines must overlap with Webster and not have too many NA's
cell_line_threshold <- 0.2
cell_line_index <- apply(prism_tmp, 1, function(x) sum(is.na(x))/length(x) < cell_line_threshold)

#Compounds not have too many missing values
compound_threshold <- 0.5
compound_index <- apply(prism_tmp[cell_line_index,],2, function(x) sum(is.na(x))/length(x) < compound_threshold)


#Impute remaining data. Takes a few mins
prism_secondary_impute <- impute_prism_na(prism_tmp[cell_line_index,compound_index])

ProjectTemplate::cache("prism_impute")


# Prepare PRISM matrix with CRISPR dictionary ------------------------------------
#remove suspension cell lines- these exhibit batch effects.


common_cl <- intersect(rownames(prism_impute), rownames(get_cell_mat(webster_depmap))) %>% 
  intersect(sample_info %>% filter(culture_type == "Adherent") %>% pull(DepMap_ID))

prism_omp_input <- prism_impute[common_cl,]

dict_input <- get_cell_mat(webster_depmap)[common_cl,] %>% 
  set_colnames(paste("V",1:webster_depmap$rank, sep = "")) %>% 
  magrittr::extract(,-c(206,215))


# Prepare orthogonal matching pursuit -------------------------------------

dict_path <- file.path(out_path, "dict_for_prism.csv")
signal_path <- file.path(out_path, "prism_signal.csv")
omp_file <- file.path(out_path, "prism_omp_output.csv")

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

cat(sys_command, sep = '\n',file=file.path(".", "cluster_scripts", paste("prism_omp","txt", sep = ".")))

#Run the script. You can also do this in the terminal
system("bash ./cluster_scripts/prism_omp.txt")

# Results -----------------------------------------------------------------
prism_omp <- read_csv(omp_file, col_names = F) %>% 
  as.matrix() %>% 
  t() %>% 
  set_rownames(colnames(prism_omp_input)) %>% 
  set_colnames(colnames(dict_input))

recon_prism <- tcrossprod(dict_input, prism_omp)

proj_results <- tibble(column_name = rownames(prism_omp),
                       Pearson = map_dbl(column_name, function(x) cor(recon_prism[,x], prism_omp_input[,x])),
                       Median_Cell_Fitness = matrixStats::colMedians(prism_omp_input)) %>% 
  left_join(secondary_meta) %>% 
  arrange(moa, name, screen_id, dose)

# heatmap -----------------------------------------------------------------


plot_moa_loadings <- function(x) {
  selected_compounds <- proj_results %>% filter(moa == x) %>% pull(column_name)
  
  tmp_omp <- prism_omp[selected_compounds,]
  
  num_of_fns <- min(sum(fn_tally > 2), 10)
  
  
  selected_fns <- colSums(tmp_omp != 0) %>% 
    enframe() %>% 
    arrange(-value) %>% 
    slice(1:num_of_fns) %>% 
    pull(name)
  
  show_numbers <- length(selected_compounds) < 40
  
  pheatmap::pheatmap(tmp_omp[,selected_fns], show_rownames = F, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", fontsize_col = 7.5,
                     annotation_row = proj_results %>% column_to_rownames("column_name") %>% select(dose, name, screen_id, Pearson, Median_Cell_Fitness) %>% unique(),
                     color = colorRampPalette(c("#BF812D", "white", "#00AEEF"))(100), main = x, display_numbers = show_numbers, number_format = "%.0f",
                     breaks = seq(-10, 10, length.out=101), filename = file.path(out_path, paste(x, "loading_heatmap.pdf", sep = "_")))
  
}

walk(moas, plot_moa_loadings)




# Recon -------------------------------------------------------------------

pita_25 <- "BRD-K75958547-238-01-0::2.5::HTS002"

cor(prism_input_cleaned[,pita_25], recon_prism[,pita_25])
plot(prism_input_cleaned[,pita_25], recon_prism[,pita_25])

orig <- cbind(prism_input_cleaned[,pita_25], recon_prism[,pita_25]) %>% scale()
tmp_heatmap <- dict_input[,c(22, 259)] %>% pheatmap::pheatmap(color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdGy")))(100) %>% rev(), breaks = seq(-0.15, 0.15, length.out=101), show_rownames = F,
                                                              cellwidth = 10, cellheight = 2, filename = "./output/figures/pivastatin_dict.pdf")

pheatmap::pheatmap(orig[tmp_heatmap$tree_row$order,], cluster_rows = F, cluster_cols = F, show_rownames = F,color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PuOr")))(100), 
                   breaks = seq(-4, 4, length.out=101), cellwidth = 10, cellheight = 2, filename = "./output/figures/pivastatin_recon.pdf")

# Dose-dependent loadings -------------------------------------------------

plot(dict_input[,259] + 4/7*dict_input[,22], prism_input_cleaned[,"BRD-K66296774-236-11-3::0.657411::HTS002"])

prism_omp[tmp_meta %>% pull(column_name),259] %>% 
  enframe("column_name", "Loading") %>% 
  left_join(tmp_meta) %>%
  mutate(Median_Fitness = matrixStats::colMedians(prism_input_cleaned[,column_name])) %>% 
  filter(name != "simvastatin") %>% 
  pivot_longer(names_to = "Mode", values_to = "Value", c(Loading, Median_Fitness)) %>% 
  mutate(Mode = factor(Mode, levels = c("Loading", "Median_Fitness"))) %>% 
  ggplot(aes(log(dose,2), Value, color)) +
  geom_point(aes(color = Mode), size = 0.5) +
  geom_col(width = 0.3, aes(fill = Mode)) +
  facet_grid(Mode~name, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Log2 Dose (uM)") +
  scale_color_manual(values=c("#2aace2", "red")) +
  scale_fill_manual(values=c("#2aace2", "red")) +
  theme(legend.position = "NA") +
  ggsave("./output/figures/loading_med_fitness_by_dose_statin.pdf", width = 6, height =2.5, device = cairo_pdf)


a+
  ggsave("./output/figures/loading_by_dose_statin.pdf", width = 5, height = 3)




#Investigations into specific cell lines

mat <- cbind(prism_input_cleaned[,tmp_meta %>% filter(name == "fluvastatin") %>% pull(column_name)], dict_input[,c(259, 22, 37, 72, 53, 208, 268)])

v37_cl <- c("ACH-000768",
            "ACH-000943",
            "ACH-000416",
            "ACH-000974",
            "ACH-000837",
            "ACH-000401")

v72_cl <- c("ACH-000040",
            "ACH-000137",
            "ACH-000211",
            "ACH-000098",
            "ACH-000469",
            "ACH-000445")

list(V37_Dependent  = v37_cl,
     V72_Dependent  = v72_cl) %>% 
  enframe("Class", "DepMap_ID") %>% 
  unnest()

mat %>% 
  as_tibble(rownames = "DepMap_ID") %>% 
  pivot_longer(names_to = "column_name", values_to = "Compound_Dependency", starts_with("BRD")) %>% 
  pivot_longer(names_to = "Function", values_to = "Fn_Dependency", starts_with("V")) %>% 
  left_join(tmp_meta) %>% 
  # left_join(list(V37_Dependent  = v37_cl,
  #                V72_Dependent  = v72_cl) %>% 
  #             enframe("Class", "DepMap_ID") %>% 
  #             unnest()) %>% 
  ggplot(aes(Compound_Dependency, Fn_Dependency, )) +
  geom_point(aes()) +
  facet_grid(Function~dose) +
  scale_color_manual(values = c("blue", "red"), na.value="black") +
  geom_smooth(method = "lm")
