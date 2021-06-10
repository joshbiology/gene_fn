library(ProjectTemplate); load.project()


load("./cache/avana_regressed_19q4.RData")
load("./cache/avana_19q4_webster.RData")
source("./munge/webster_depmap.R")


get_gene_mat(webster_depmap)[convert_genes("TP53", "symbol", "cds_id"),] %>% enframe()

#Loadings are in units of variance.
norm_info <- tibble(Gene = colnames(avana_19q4_webster),
       Norm_Orig = map_dbl(Gene, ~norm(avana_19q4_webster[,.], type="2")),
       Norm_Recon = map_dbl(Gene, ~norm(recon_depmap[,.], type="2")),
       Ratio = Norm_Recon/Norm_Orig,
       Pearson = map_dbl(Gene, ~cor(avana_19q4_webster[,.], recon_depmap[,.])),
       Recon_Error = map_dbl(Gene, ~norm(avana_19q4_webster[,.] - recon_depmap[,.], type = "2"))) %>% 
  arrange(Ratio)


#We can correct loadings such that the median (non-zero) loading is equal to the median gene norm in the avana dataset.


orig_norms <- tibble(Gene = colnames(avana_regressed_19q4),
       Norm_Orig = map_dbl(Gene, ~norm(avana_regressed_19q4[,.], type="2")))

hist(orig_norms$Norm_Orig)

hist(orig_norms %>% filter(Gene %in% colnames(avana_19q4_webster)) %>% pull(Norm_Orig))

median(orig_norms$Norm_Orig/sqrt(675))

matrixStats::colSds(avana_regressed_19q4) %>% median()
matrixStats::colSds(avana_19q4_webster) %>% median()

matrixStats::rowSds(avana_regressed_19q4)
matrixStats::rowSds(avana_19q4_webster) %>% summary()
matrixStats::colSds(avana_19q4_webster) %>% summary()

get_cell_mat(webster_depmap)[,1] %>% norm(type = "2")
(get_cell_mat(webster_depmap)[,1]*26)%>% norm(type = "2")
(get_cell_mat(webster_depmap)[,1]*26)%>% sd()


hist(norm_info$Norm_Orig)
median(norm_info$Norm_Orig/sqrt(675))

webster_loadings <- get_gene_mat(webster_depmap) %>% abs() %>% 
  as_tibble(rownames = "Gene") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V")) %>% 
  filter(abs(Loading) > 0) %>% 
  group_by(Gene) %>% 
  arrange(Gene, desc(abs(Loading))) %>% 
  mutate(Rank = 1:4) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Rank, values_from = c("Function", "Loading"))


#https://pballew.blogspot.com/2010/07/standard-deviation-as-distance.html
webster_loadings_stdev_unit <- (get_gene_mat(webster_depmap)/sqrt(675)) %>% 
  abs %>% 
  as_tibble(rownames = "Gene") %>% 
  pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V")) %>% 
  filter(abs(Loading) > 0) %>% 
  group_by(Gene) %>% 
  arrange(Gene, desc(abs(Loading))) %>% 
  mutate(Rank = 1:4) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Rank, values_from = c("Function", "Loading"))

webster_loadings %>% 
  left_join(norm_info) %>% View()

webster_loadings_stdev_unit %>% 
  left_join(norm_info) %>% View()

median(get_gene_mat(webster_depmap) %>% abs() %>% 
         as_tibble(rownames = "Gene") %>% 
         pivot_longer(names_to = "Function", values_to = "Loading", starts_with("V")) %>% 
         filter(abs(Loading) > 0) %>% pull(Loading))

recon_depmap[,convert_genes("TP53", "symbol", "cds_id")] %>% norm(type = "2")

avana_regressed_19q4 %>% std()

colMeans(avana_regressed_19q4) %>% summary()
matrixStats::colVars(avana_regressed_19q4) %>% summary()
tmp <- var(avana_regressed_19q4 %>% as.numeric())


avana_19q4_webster %>% aaply(2,norm_vec) %>% summary()


get_cell_mat(webster_depmap) %>% 
  colMeans() %>% summary()


get_cell_mat(webster_depmap) %>% 
  matrixStats::colVars() %>% summary()

norm_vec <- function(x) sqrt(sum(x^2))

plot(avana_19q4_webster[,1],recon_depmap[,1] * 1.01/.2848)

plot(avana_19q4_webster[,2921],recon_depmap[,2921] )

plot(avana_19q4_webster[,2921],get_cell_mat(webster_depmap)[,38] * 60)

plot(avana_19q4_webster[,2],recon_depmap[,2] *  1.017333/.115983)


tmp1 <- avana_19q4_webster %>% aaply(2, function(x)sqrt(var(x)))
tmp2 <- recon_depmap %>% aaply(2, function(x)sqrt(var(x)))
tmp3 <- get_cell_mat(webster_depmap)  %>% aaply(2, function(x)sqrt(var(x)))


((avana_19q4_webster %>% aaply(2, norm_vec))/(recon_depmap %>% aaply(2, norm_vec)) ) %>% summary()

((avana_19q4_webster %>% aaply(2, function(x)sqrt(var(x))))/(recon_depmap %>% aaply(2, function(x)sqrt(var(x)))) ) %>% summary()
