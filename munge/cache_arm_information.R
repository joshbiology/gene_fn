#cache arm information

library(tidyverse)
library(magrittr)
library(biomaRt)

#ensembl <- useMart("ensembl")
#ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl, mirror="useast")

#listAttributes(ensembl)

new_config <- httr::config(ssl_verifypeer = FALSE)
httr::set_config(new_config, override = FALSE)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

attr_tmp <- c("entrezgene_id", "hgnc_symbol", "chromosome_name", "band")

chrom_names <- c(1:22, "X", "Y")

arm_df <- getBM(attr_tmp, values = 'entrezgene_id', mart = ensembl) %>% 
  dplyr::filter(chromosome_name %in% chrom_names, !is.na(entrezgene_id)) %>% 
  dplyr::mutate(arm = str_sub(band, 1, 1)) %>% 
  dplyr::rename(entrezgene = entrezgene_id, symbol = hgnc_symbol) %>% 
  dplyr::mutate(entrezgene = as.character(entrezgene)) %>%
  as_tibble()

ProjectTemplate::cache("arm_df")
