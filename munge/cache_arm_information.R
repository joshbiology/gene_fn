#cache arm information

library(tidyverse)
library(magrittr)
library(biomaRt)

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

listAttributes(ensembl)

attr_tmp <- c("entrezgene_id", "hgnc_symbol", "chromosome_name", "band")

chrom_names <- c(1:22, "X", "Y")

arm_df <- getBM(attr_tmp, values = 'entrezgene_id', mart = ensembl) %>% 
  filter(chromosome_name %in% chrom_names, !is.na(entrezgene_id)) %>% 
  mutate(arm = str_sub(band, 1, 1)) %>% 
  rename(entrezgene = entrezgene_id, symbol = hgnc_symbol) %>% 
  mutate(entrezgene = as.character(entrezgene))
  as_tibble()

ProjectTemplate::cache("arm_df")