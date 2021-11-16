library(ProjectTemplate); load.project()

require(future)
require(purrr)
plan(multisession)


out_path <- file.path(".", "output", "04-depmap_postprocess")
create_output_folder(out_path)

# Functions ---------------------------------------------------------------



graph_objective_precomp <- function(factorized_mat, L_mat) {
  sum(diag(Rfast::Crossprod(factorized_mat,L_mat) %>%
             Rfast::mat.mult(factorized_mat)))
}


import_graphDL  <- function(x, input = avana_19q4_webster) {
  tmp <- R.matlab::readMat(x)

  graphdl_to_factorized(import_graphdl(tmp), colnames(input), rownames(input))
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
    dplyr::mutate(F_Norm = f_norms,
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
  dplyr::mutate(Gene_Lap = log(Gene_Lap, 10)) %>%
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>%
  dplyr::filter(K < 600) %>%
  ggplot(aes(x = K, y = Value, group = T_param, color = factor(T_param))) +
  geom_point() +
  geom_line() +
  facet_wrap(~Metric, scales ="free_y", ncol = 2)
  ggsave(file.path(out_path,"depmap_grid_metrics.pdf"), width = 7, height = 3, device = cairo_pdf)


# Wide ------------------------------------------------------------------

#Wide factorization output: many values of K, one value of T, one random seed
mat_paths <- list("./data/interim/matlab/depmap_wide") %>%
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

depmap_wide_output <- map(mat_paths, import_graphDL)

tictoc::tic()
depmap_wide_metrics_df <- generate_metrics(depmap_wide_output)
tictoc::toc()

write_tsv(depmap_wide_metrics_df, path = file.path(out_path, "depmap_wide_metrics.tsv"))
#If loading from saved TSV:
#depmap_wide_metrics_df <- read_tsv(file.path(out_path, "depmap_wide_metrics.tsv"))


depmap_wide_metrics_df %>%
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>%
  ggplot(aes(x = K, y = Value)) +
  geom_jitter() +
  geom_smooth() +
  geom_vline(xintercept = 220) +
  xlim(c(0, 600)) +
  facet_grid( Metric ~Num_Neigh_Gene, scales ="free_y")
  ggsave(file.path(out_path,"depmap_wide_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)

depmap_wide_metrics_df %>%
  arrange(K) %>%
  dplyr::mutate(Diff_F = c(0,diff(F_Norm)),
         Diff_Gene_L = c(0,diff(Gene_Lap))) %>%
  slice(-1) %>%
  pivot_longer(names_to = "Metric", values_to = "Value", c("Diff_F", "Diff_Gene_L")) %>%
  dplyr::filter(K < 500) %>%
  ggplot(aes(x = K, y = Value))+
  geom_vline(xintercept = 220) +
  xlim(c(0, 500)) +
  geom_jitter() +
  geom_smooth() +
  facet_wrap( Metric ~Num_Neigh_Gene, scales ="free_y")
  ggsave(file.path(out_path,"depmap_wide_marginal_metrics_T=4.pdf"), width = 7, height = 3.5, device = cairo_pdf)

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

#For reproducibility:
#depmap_deep_metrics_df <- read_tsv(file.path(out_path, "depmap_deep_metrics.tsv"))
depmap_deep_metrics_df %>%
  pivot_longer(names_to = "Metric", values_to = "Value", c(F_Norm, Gene_Lap)) %>%
  ggplot(aes(x = K, y = Value)) +
  geom_jitter() +
  geom_smooth() +
  facet_grid( Metric ~Num_Neigh_Gene, scales ="free_y")
  ggsave(file.path(out_path,"depmap_deep_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)

grouped_metrics <- depmap_deep_metrics_df %>%
  group_by(K) %>% summarize(Mean_F_Norm = mean(F_Norm), Mean_Gene_Lap = mean(Gene_Lap))

grouped_metrics %>%
  dplyr::mutate(Diff_F = c(0,diff(Mean_F_Norm)),
         Diff_Gene_L = c(0,diff(Mean_Gene_Lap))) %>%
  slice(-1) %>%
  pivot_longer(names_to = "Metric", values_to = "Value", c("Diff_F", "Diff_Gene_L")) %>%
  dplyr::filter(K < 500) %>%
  ggplot(aes(x = K, y = Value))+
  geom_jitter() +
  geom_smooth() +
  facet_wrap(  ~Metric, ncol = 1, scales ="free_y")
  ggsave(file.path(out_path,"depmap_deep_marginal_metrics_T=4.pdf"), width = 6, height = 10, device = cairo_pdf)


#Consistency between dictionary atom runs.
random_seed_crosscor <- cor(depmap_deep_output[[71]]$cell_mat, depmap_deep_output[[72]]$cell_mat)

init_final_dict_crosscor <- cor(avana_19q4_webster[, attr(depmap_deep_output[[71]], "extras")$medoids[,1]], depmap_deep_output[[71]]$cell_mat)

pheatmap::pheatmap(tmp, cluster_cols = F, cluster_rows = F, breaks = seq(-1,1, length.out = 100))
pheatmap::pheatmap(tmp2, cluster_cols = F, cluster_rows = F, breaks = seq(-1, 1, length.out = 100), show_rownames = F)

#Print medoids names
colnames(avana_19q4_webster)[attr(depmap_deep_output[[71]], "extras")$medoids[,1]]

tibble(Init_Final = init_final_dict_crosscor %>% diag(),
       Random_Seed = random_seed_crosscor%>% diag()) %>%
  pivot_longer(names_to = "Type", values_to = "Cross_Cor", everything()) %>%
  ggplot(aes(Cross_Cor)) +
  geom_histogram() +
  facet_wrap(~Type)
  ggsave(file.path(out_path, "consistency_hist.pdf"), width = 4, height = 1.5)


median(random_seed_crosscor %>% diag)
median(init_final_dict_crosscor %>% diag)

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



# Section II: Denoising  --------------------------------------------------
#This is not really a post-processing step but I'm putting it here next to the transferability analysis, which is.

#Create test train split.
all_genes <- colnames(avana_19q4_webster)

test_genes <- sample(all_genes, length(all_genes)/4)
training_genes <- setdiff(all_genes, test_genes)

train_test_df <- list(Test = test_genes, Train = training_genes) %>%
  enframe("Group","CDS_ID") %>%
  unnest(CDS_ID)

train_test_df %>%
  write_tsv(file.path(out_path, paste(Sys.Date(), "train_test_genes.csv", sep  = "_")))

#For reproducibility:
#train_test_df <- read_tsv(file.path(out_path, "2021-07-01_train_test_genes.csv"))
#test_genes <- train_test_df %>% filter(Group == "Test") %>% pull(CDS_ID)
#training_genes <-  train_test_df %>% filter(Group == "Train") %>% pull(CDS_ID)


#Create train and test matrices.
avana_19q4_webster_test <- avana_19q4_webster[,filter(train_test_df, Group == "Test")$CDS_ID]
avana_19q4_webster_train  <- avana_19q4_webster[,filter(train_test_df, Group == "Train")$CDS_ID]


# #Create noise matrices.
# noise_25 <-  rnorm(length(avana_19q4_webster), sd = 0.25) %>% matrix(nrow = dim(avana_19q4_webster)[1]) %>% set_colnames(colnames(avana_19q4_webster))
# noise_50 <-  rnorm(length(avana_19q4_webster), sd = 0.50) %>% matrix(nrow = dim(avana_19q4_webster)[1]) %>% set_colnames(colnames(avana_19q4_webster))
# noise_100 <-  rnorm(length(avana_19q4_webster), sd = 1) %>% matrix(nrow = dim(avana_19q4_webster)[1]) %>% set_colnames(colnames(avana_19q4_webster))
# noise_150 <-  rnorm(length(avana_19q4_webster), sd = 1.5) %>% matrix(nrow = dim(avana_19q4_webster)[1]) %>% set_colnames(colnames(avana_19q4_webster))
#
# #Exporting noise matrices for reproducibility.
# export_for_matlab(noise_25, file.path(out_path, "noise_25.csv"))
# export_for_matlab(noise_50, file.path(out_path, "noise_50.csv"))
# export_for_matlab(noise_100, file.path(out_path, "noise_100.csv"))
# export_for_matlab(noise_150, file.path(out_path, "noise_150.csv"))

#Import for reproducibility:
noise_25 <- read_csv(file.path(out_path, "noise_25.csv"), col_names =  colnames(avana_19q4_webster)) %>% as.matrix()
noise_50 <- read_csv(file.path(out_path, "noise_50.csv"), col_names =  colnames(avana_19q4_webster)) %>% as.matrix()
noise_100 <- read_csv(file.path(out_path, "noise_100.csv"), col_names =  colnames(avana_19q4_webster)) %>% as.matrix()
noise_150 <- read_csv(file.path(out_path, "noise_150.csv"), col_names =  colnames(avana_19q4_webster)) %>% as.matrix()

#Add noise.
noisy_input_25 <- avana_19q4_webster+noise_25
noisy_input_50 <- avana_19q4_webster+noise_50
noisy_input_100 <- avana_19q4_webster+noise_100
noisy_input_150 <- avana_19q4_webster+noise_150

plot(avana_19q4_webster[,1], noisy_input_25[,1])
plot(avana_19q4_webster[,1], noisy_input_50[,1])
plot(avana_19q4_webster[,1], noisy_input_100[,1])
plot(avana_19q4_webster[,1], noisy_input_150[,1])

#Export noisy data.
#Test
export_for_matlab(avana_19q4_webster[,test_genes]+noise_25[,test_genes], file.path(out_path, "avana_19q4_webster_test_noise_25.csv"))
export_for_matlab(avana_19q4_webster[,test_genes]+noise_50[,test_genes], file.path(out_path, "avana_19q4_webster_test_noise_50.csv"))
export_for_matlab(avana_19q4_webster[,test_genes]+noise_100[,test_genes], file.path(out_path, "avana_19q4_webster_test_noise_100.csv"))
export_for_matlab(avana_19q4_webster[,test_genes]+noise_150[,test_genes], file.path(out_path, "avana_19q4_webster_test_noise_150.csv"))

#Train
export_for_matlab(avana_19q4_webster[,training_genes]+noise_25[,training_genes], file.path(out_path, "avana_19q4_webster_train_noise_25.csv"))
export_for_matlab(avana_19q4_webster[,training_genes]+noise_50[,training_genes], file.path(out_path, "avana_19q4_webster_train_noise_50.csv"))
export_for_matlab(avana_19q4_webster[,training_genes]+noise_100[,training_genes], file.path(out_path, "avana_19q4_webster_train_noise_100.csv"))
export_for_matlab(avana_19q4_webster[,training_genes]+noise_150[,training_genes], file.path(out_path, "avana_19q4_webster_train_noise_150.csv"))

#Export original matrices.
export_for_matlab(avana_19q4_webster_test, file.path(out_path, "avana_19q4_webster_test_noise_0.csv"))
export_for_matlab(avana_19q4_webster_train, file.path(out_path, "avana_19q4_webster_train_noise_0.csv"))

#Write bash script (This step is performed in 00-webster_scripts.R)

#Execute batch run locally
system("bash ./cluster_scripts/denoise.txt")



#Import datasets
mat_paths <- list("./data/interim/matlab/depmap_denoise") %>%
  map(~list.files(., pattern = "*.mat", full.names = T)) %>% unlist()

depmap_denoise_output <- map(mat_paths, ~import_graphDL(., avana_19q4_webster_train))

denoise_params_df <- extract_params(depmap_denoise_output)

#Write dictionaries to file
walk(1:length(depmap_denoise_output), function(x) {
  dict <- depmap_denoise_output[[x]] %>% get_cell_mat %>%
    export_for_matlab(file.path(out_path, sprintf("%s_%s.csv", list.files("./data/interim/matlab/depmap_denoise", pattern = "*.mat", full.names = F)[x], x)))
})

#write script to perform OMP using noisy test data and the noisy trained dictionary.
omp_params <- tibble(dict_path = map_chr(1:length(depmap_denoise_output), ~file.path(out_path, sprintf("%s_%s.csv", list.files("./data/interim/matlab/depmap_denoise", pattern = "*.mat", full.names = F)[.], .))),
       signal_path = c(rep(file.path(out_path, "avana_19q4_webster_test_noise_0.csv"), 5),
                       rep(file.path(out_path, "avana_19q4_webster_test_noise_100.csv"), 5),
                       rep(file.path(out_path, "avana_19q4_webster_test_noise_150.csv"), 5),
                       rep(file.path(out_path, "avana_19q4_webster_test_noise_25.csv"), 5),
                       rep(file.path(out_path, "avana_19q4_webster_test_noise_50.csv"), 5)),
       omp_file = file.path(out_path, outer(paste(c(0, 100, 150, 25, 50), "denoise_omp", sep = "_"), 1:5, paste, sep = "_") %>% t() %>%  as.vector() %>% paste(".csv", sep = "")))

sys_command <- map_chr(1:nrow(omp_params), function(x){

  sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
          omp_params[[x,1]],
          omp_params[[x,2]],
          omp_params[[x,3]],
          4)
})

cat(sys_command, sep = '\n',file=file.path(".", "cluster_scripts", paste("denoise_omp","txt", sep = ".")))

#Execute
system("bash ./cluster_scripts/denoise_omp.txt")

#Read OMP matrices
denoise_omp <- list.files(out_path, pattern = "*denoise_omp*", full.names = T) %>%
  map(~read_csv(., col_names = test_genes) %>% as.matrix)

test_recon <- map2(depmap_denoise_output, denoise_omp, ~ get_cell_mat(.x) %*% .y)

noisy_pearson_df <- map_dfr(test_recon, function(x) tibble(Gene = test_genes,
                                                       Pearson_Recon = map_dbl(1:length(test_genes), ~cor(avana_19q4_webster_test[,.], x[,.]))), .id = "ID") %>%
  left_join(tibble(ID = as.character(1:25), Noise_Level = rep(c(0, 100, 150, 25, 50), each = 5), Seed = rep(1:5, 5))) %>%
  group_by(Noise_Level, Gene) %>%
  summarize(Mean_Pearson = mean(Pearson_Recon))


noisy_pearson_df %>%
  group_by(Noise_Level) %>%
  dplyr::mutate(Mean = mean(Mean_Pearson)) %>%
  ungroup() %>%
  dplyr::mutate(Noise_Level = factor(Noise_Level, labels = c("σ=0", "σ=0.25", "σ=0.5", "σ=1", "σ=1.5"))) %>%
  ggplot(aes(x = Mean_Pearson)) +
  geom_density(fill = "gray90") +
  geom_vline(aes(xintercept = Mean), color = "red", linetype = "dashed") +
  facet_wrap(~Noise_Level, ncol = 1)
  ggsave(file.path(out_path, "denoise.pdf"), width =3, height = 4, device = cairo_pdf)



# Section III: Transferability  --------------------------------------------------

#use dense factorized output to generate Sanger recoveries.
source("./munge/webster_depmap.R")

load("./cache/sanger_regressed_19q4.RData")

common_genes <- intersect(colnames(sanger_regressed_19q4), get_gene_mat(webster_depmap) %>% rownames)

common_cl <- intersect(rownames(sanger_regressed_19q4), get_cell_mat(webster_depmap) %>% rownames)

shuffle_1 <- sample(common_cl)
shuffle_2 <-  sample(common_cl)

common_dict <- get_cell_mat(webster_depmap)[common_cl,]
sanger_signal <-  sanger_regressed_19q4[common_cl, common_genes]

#Export
export_for_matlab(common_dict, file.path(out_path, "sanger_dict.csv"))
export_for_matlab(common_dict[shuffle_1,], file.path(out_path, "sanger_dict_shuffled.csv"))
export_for_matlab(sanger_signal, file.path(out_path, "sanger_signal.csv"))
export_for_matlab(avana_19q4_webster[common_cl, common_genes], file.path(out_path, "avana_signal.csv"))


sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
        file.path(out_path, "sanger_dict.csv"),
        file.path(out_path, "sanger_signal.csv"),
        file.path(out_path, "sanger_omp.csv"),
        4) %>%
  cat(file = file.path(".", "cluster_scripts", "sanger.txt"))

system("bash ./cluster_scripts/sanger.txt")

sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
        file.path(out_path, "sanger_dict_shuffled.csv"),
        file.path(out_path, "sanger_signal.csv"),
        file.path(out_path, "sanger_omp_shuffled.csv"),
        4) %>%
  cat(file = file.path(".", "cluster_scripts", "sanger_shuffled.txt"))

system("bash ./cluster_scripts/sanger_shuffled.txt")

sprintf("/Applications/MATLAB_R2019b.app/bin/matlab -batch  \"OMP_wrapper('%s', '%s', '%s', %d);exit\"",
        file.path(out_path, "sanger_dict.csv"),
        file.path(out_path, "avana_signal.csv"),
        file.path(out_path, "avana_omp.csv"),
        4) %>%
  cat(file = file.path(".", "cluster_scripts", "avana.txt"))

system("bash ./cluster_scripts/avana.txt")

avana_omp <- read_csv(file.path(out_path, "avana_omp.csv"), col_names = common_genes) %>% as.matrix
sanger_omp <- read_csv(file.path(out_path, "sanger_omp.csv"), col_names = common_genes) %>% as.matrix
sanger_omp_shuffled <- read_csv(file.path(out_path, "sanger_omp_shuffled.csv"), col_names = common_genes) %>% as.matrix

recon_avana <- common_dict%*% avana_omp

recon_sanger <- common_dict%*% sanger_omp

recon_sanger_shuffled <- common_dict[shuffle_1,]%*% sanger_omp_shuffled

recon_df <- tibble(Gene= common_genes,
                   Sanger_Recon_Pearson = map_dbl(common_genes, ~cor(sanger_signal[,.], recon_sanger[,.])),
                   Shuffled_Recon_Pearson = map_dbl(common_genes, ~cor(sanger_signal[common_cl,.], recon_sanger_shuffled[common_cl,.])),
                   Avana_Recon_Pearson =  map_dbl(common_genes, ~cor(avana_19q4_webster[common_cl,.], recon_avana[common_cl,.])))

recon_df %>%
  pivot_longer(names_to = "Metric", values_to = "Value", -Gene) %>%
  group_by(Metric) %>%
  dplyr::mutate(Mean = mean(Value)) %>%
  ungroup() %>%
  ggplot(aes(Value)) +
  geom_density(fill = "gray90") +
  geom_vline(aes(xintercept = Mean), color = "red", linetype = "dashed") +
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  ggsave(file.path(out_path, "density_recon_sanger.pdf"), width =3, height = 3)
