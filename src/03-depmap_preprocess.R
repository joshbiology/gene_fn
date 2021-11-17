#import 19q4 data -- public and ready to use
#output the gene selection features for post-model prediction.


#19q4 public
#but the fitness genes I use are done as in my original filtering.
#Notes for Avana
#1. Arm correction - none
#2. MON artifact - 8 worst cell lines were removed. Other cell lines have those genes as NA
#3. pDNA batch effect - taken care of
#4. Multiguide correction - CERES uses all guide data, but genes whose guides are multitargeting are removed entirely.
#5. Genes without X chrom copy number data - NA'ed

#Notes for Sanger
#1. Same as 5/2019 release by Sanger

#Protocol
#1. Remove the 12 MON cell lines.
#2. Remove the remaining X-chrom genes that are missing data in 47 cell lines.
#3. Perform arm correction and store metrics.
#4. Scale and regress out NNMD, store metrics.


library(ProjectTemplate); load.project()

out_path <- file.path(".", "output", "03-depmap_preprocess")
create_output_folder(out_path)



# Functions -------------------------------------------------------
#Gene effect matrix - rows = cell lines, cols = genes. Gene names = CDS_ID, format %HUGO (%Entrez)

query_na <- function(gene_effect_matrix) {
  genes <- matrixStats::colSums2(is.na(gene_effect_matrix)) %>%
    set_names(colnames(gene_effect_matrix)) %>%
    enframe() %>%
    count(value) %>%
    set_colnames(c("Number of NA's", "Number of genes"))


  cell_lines <- matrixStats::rowSums2(is.na(gene_effect_matrix)) %>%
    set_names(rownames(gene_effect_matrix)) %>%
    enframe() %>%
    count(value) %>%
    set_colnames(c("Number of NA's", "Number of cell lines"))

  return(list(Genes = genes, Cell_Lines = cell_lines))
}


remove_nuisance_genes <- function(gene_effect_matrix) {
  #These are manually defined nusiance gene sets from HUGO. Ex: olf receptors (Boyle et al. 2018)

  nusiance_gene_df <- list.files("./data/raw/nusiance_genesets", pattern = "tsv", full.names = T) %>%
    map_df(read_tsv)

  gene_effect_matrix[, !(convert_cds_to_entrez(colnames(gene_effect_matrix)) %in% nusiance_gene_df$`NCBI Gene ID`)]

}

center_cell_lines <- function(gene_effect_matrix) {
  means <- rowMeans(gene_effect_matrix)
  sweep(gene_effect_matrix, 1, means)
  #apply(gene_effect_matrix, 1, function(y) y - mean(y)) %>% t()
}

regress_nnmd <- function(gene_effect_matrix, cell_line_names, nnmd_values) {
  nnmd <- nnmd_values
  names(nnmd) <- cell_line_names


  plyr::aaply(gene_effect_matrix, 2, function(x){
    lm(x~ nnmd[rownames(gene_effect_matrix)])$residuals
  }) %>% scale(center = F, scale = T) %>% t()
}

# Import and remove NA's from 19q4 -------------------------------------------------------------
gene_effect_19q4 <- read_tsv("./data/raw/depmap/public-19q4_v23-achilles-gene-effect.tsv") %>%
  column_to_rownames("...1") %>%
  as.matrix()

query_na(gene_effect_19q4)

#The cell lines with > 1200 NA's represent bad batches of cell lines. Remove them
bad_batch_threshold <- 1200
tmp <- gene_effect_19q4[matrixStats::rowSums2(is.na(gene_effect_19q4)) < bad_batch_threshold,]
query_na(tmp)

#The remaining genes with NA's are sex chrom. genes with missing copy number data. Remove them
tmp2 <-tmp %>%
  remove.cols.any.nas()

query_na(tmp2)

#Remove nusiance genes
avana_tmp <- remove_nuisance_genes(tmp2)


# Import and remove NA's from Sanger data ------------------------------------------------------

#Import and remove any genes with NA's
sanger <- read_tsv("./data/raw/depmap/sanger-crispr-project-score-_v4-gene-effect.tsv") %>%
  column_to_rownames("...1") %>%
  as.matrix() %>%
  remove.cols.any.nas()

query_na(sanger)

#Remove nusiance genes
sanger_tmp <- remove_nuisance_genes(sanger)


# Arm correct datasets -----------------------------------------------------

avana_arm_corrected <- arm_correct(avana_tmp)
sanger_arm_corrected <- arm_correct(sanger_tmp)


# Center cell lines  -------------------------------------------------------
avana_centered <- center_cell_lines(avana_arm_corrected)
sanger_centered <- center_cell_lines(sanger_arm_corrected)

# Regress out NNMD and scale cell lines ------------------------------------
#For Avana
sample_info_19q4 <- read_tsv("./data/raw/depmap/public-19q4_v23-sample-info.tsv")
avana_regressed <- regress_nnmd(avana_centered, sample_info_19q4$DepMap_ID, sample_info_19q4$cell_line_NNMD)


#For Sanger
nnmd_sanger <- read_csv("./data/raw/depmap/score_v4_nnmd.csv", col_names = c("DepMap_ID", "cell_line_NNMD"))
sanger_regressed <- regress_nnmd(sanger_centered, nnmd_sanger$DepMap_ID, nnmd_sanger$cell_line_NNMD)


# Top genes by variance -----------------------------------------------
nonessentials <- read_tsv("./data/raw/depmap/public-19q4_v23-nonessentials.tsv")$gene
common.essentials <- read_tsv("./data/raw/depmap/public-19q4_v23-common-essentials.tsv")$gene

fitness_list <- list(Avana = avana_regressed,
                     Sanger = sanger_regressed)


gene_vars <- map_df(fitness_list, function(mat) mat %>%
                      matrixStats::colVars(na.rm = T) %>%
                      set_names(colnames(mat)) %>%
                      enframe("Gene", "Var") %>%
                      dplyr::mutate(Rank = order(order(Var))), .id = "Dataset")

gene_vars %>%
  ggplot(aes(Var)) +
  geom_histogram(bins = 100) +
  geom_freqpoly(bins = 100) +
  scale_x_log10() +
  ggtitle("Distribution of gene variances (log10)") +
  theme_minimal() +
  facet_wrap(~Dataset, ncol = 1)
  ggsave(file.path(out_path,"log_gene_var_histogram_19q4.pdf"), width = 4, height = 3)

gene_vars %>%
  dplyr::select(Gene, Dataset, Var) %>%
  pivot_wider(names_from = "Dataset", values_from = "Var") %>%
  ggplot(aes(Avana, Sanger)) +
  geom_point(alpha = 0.2, size = 0.5) +
  scale_x_log10(limits = c(0.1, NA)) +
  scale_y_log10() +
  geom_density_2d(binwidth = 0.2, alpha = 0.75, color = "cyan") +
  ggtitle("Distribution of gene variances (log10)") +
  theme_minimal()
  ggsave(file.path(out_path,"gene_variance_scatter_19q4.png"), width = 4.5, height = 4)

gene_vars %>%
  ggplot(aes(sample = Var)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~Dataset) +
  ggtitle("QQplot of gene variances")

#Take genes that deviate from the qqplot
#https://mgimond.github.io/ES218/Week06a.html
qq_data <- gene_vars %>%
  group_by(Dataset) %>%
  arrange(Dataset, Rank) %>%
  dplyr::mutate(Quantile = ppoints(length(Rank))) %>%
  dplyr::mutate(Theoretical = qnorm(Quantile)) %>%
  ungroup()

qq_fit_line_avana <- lm(Var ~ Theoretical, data = qq_data %>%
                          dplyr::filter(Theoretical < 0, Theoretical > -3, Dataset == "Avana"))

qq_fit_line_sanger <- lm(Var ~ Theoretical, data = qq_data %>%
                           dplyr::filter(Theoretical < 0, Theoretical > -3, Dataset == "Sanger"))

#Rank cutoff determines which genes are included in the linear regression
min_rank <- 5000
#min_diff determines how far variance should be from the line to be counted as "high" variance
min_diff <- 0.3

qq_data <- qq_data %>%
  dplyr::mutate(Predicted = c(predict(qq_fit_line_avana,
                               split(qq_data, qq_data$Dataset)$Avana %>% dplyr::select(Theoretical)),
                       predict(qq_fit_line_sanger,
                               split(qq_data, qq_data$Dataset)$Sanger %>% dplyr::select(Theoretical))),
         Diff = Var - Predicted,
         High_Var = Diff > min_diff & Rank > min_rank)

qq_data %>%
  ggplot(aes(Theoretical, Var)) +
  geom_point(size = 0.5) +
  geom_abline(slope = qq_fit_line_avana$coefficients[2], intercept = qq_fit_line_avana$coefficients[1],
              linetype = "dashed") +
  ggtitle("QQplot of Gene Variances") +
  scale_color_manual(values = c("FALSE" = 'gray', "TRUE" = 'red')) +
  theme_minimal() +
  facet_wrap(~Dataset)
  ggsave(file.path(out_path,"gene_qqplot_19q4.pdf"), width = 6, height = 3)

# Gene confidence ---------------------------------------------------------
#http://archive.today/2021.03.22-122633/https://cancerdatascience.org/blog/posts/gene_confidence_blog/

gene.confidence.features <- read_tsv("./data/raw/depmap/gene-confidence-internal-20q4_v2-gene-confidence-features.tsv") %>%
  dplyr::select(-1) %>%
  rename(Gene = gene)

qq_data %>% dplyr::mutate(entrezgene = convert_cds_to_entrez(Gene)) %>%
  left_join(gene.confidence.features %>%
              dplyr::mutate(entrezgene = convert_cds_to_entrez(Gene)) %>%
              dplyr::select(-Gene), by = "entrezgene") %>%
  ggplot(aes(Var, confidence)) +
  geom_point(alpha = 0.25) +
  scale_x_log10() +
  facet_wrap(~Dataset)

gene_selection_features <- qq_data %>%
  dplyr::mutate(entrezgene = convert_cds_to_entrez(Gene)) %>%
  left_join(gene.confidence.features %>%
              dplyr::mutate(entrezgene = convert_cds_to_entrez(Gene)) %>%
              dplyr::select(-Gene), by = "entrezgene")

# Cache -------------------------------------------------------------------
avana_regressed_19q4 <- avana_regressed
sanger_regressed_19q4 <- sanger_regressed
gene_selection_features_19q4 <- gene_selection_features
ProjectTemplate::cache("avana_regressed_19q4")
ProjectTemplate::cache("sanger_regressed_19q4")
ProjectTemplate::cache("gene_selection_features_19q4")


# Feature selection: max cor ---------------------------------------------
#List of genes with above-background variance
total_genes <- gene_selection_features_19q4 %>%
  dplyr::filter(Dataset=="Avana", Diff > 0.05)

cor_tmp <- coop::pcor(avana_regressed_19q4[,total_genes$Gene])

diag(cor_tmp) <- NA

total_genes <- total_genes %>%
  dplyr::mutate(Max_Cor = matrixStats::rowMaxs(abs(cor_tmp), na.rm = T))

write_csv(total_genes, file.path(out_path,"genes_above_background_var_features_19q4.csv"))



# Feature selection: expression addition essentials -----------------------
#These genes are essential genes that exhibit expression additions to a differentially expressed paralogue, and therefore have low fitness correlations
#with other genes due to the specificity of their dependency. Idenfied manually, as I couldn't find a good cutoff for them.
high_var_outliers <- c("BCL2L1 (598)", "NXT1 (29107)", "SNRPB2 (6629)", "FBXW11 (23291)", "SLC2A1 (6513)")

# Feature selection: thresholding on variance, max_cor, conf --------------
correlation_threshold <- 0.275
diff_threshhold <- 0.25
conf_threshold <- 0.5


g1 <- total_genes %>%
  dplyr::filter(confidence > conf_threshold) %>%
  ggplot(aes(Diff, Max_Cor, color = confidence))+ scale_x_log10() + geom_point(alpha = 0.5) +
  geom_hline(yintercept = correlation_threshold) +
  geom_vline(xintercept = diff_threshhold) +
  scale_color_viridis_c(limits = c(0,1))


g2 <- total_genes %>%
  dplyr::filter(confidence <= conf_threshold) %>%
  ggplot(aes(Diff, Max_Cor, color = confidence))+ scale_x_log10() + geom_point(alpha = 0.5) +
  geom_hline(yintercept = correlation_threshold) +
  geom_vline(xintercept = diff_threshhold) +
  scale_color_viridis_c(limits = c(0,1))


require(ggExtra)
g3 <- ggMarginal(g1,type="histogram")

ggsave(g3, filename = file.path(out_path, "high_conf_perturbation.pdf"), width = 6.5, height = 6, device = cairo_pdf)
g4 <- ggMarginal(g2,type="histogram")
ggsave(g4, filename = file.path(out_path, "low_conf_perturbation.pdf"), width = 6.5, height = 6, device = cairo_pdf)



# Filtered Webster input -------------------------------------------------------------------
avana_19q4_webster <- avana_regressed_19q4[,total_genes %>%
                                                       dplyr::filter(confidence > conf_threshold,
                                                              Max_Cor > correlation_threshold,
                                                              Diff > diff_threshhold) %>%
                                             pull(Gene) %>%
                                             setdiff(high_var_outliers)]


ProjectTemplate::cache("avana_19q4_webster")

#Export for matlab -------------------------------------------------------------------
export_for_matlab(avana_19q4_webster, out_file = "./output/03-depmap_preprocess/avana_19q4_webster.csv")
