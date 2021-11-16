#cache arm information
#Goal: To query chrom arm information from biomaRt. 
#The output of this file is stored as a data file, arm_df.tsv, that is part of the data repository on FigShare.
#So this script is never run for new users, but is maintained here for reproducibility.

library(tidyverse)
library(magrittr)
library(biomaRt)

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

write_tsv(arm_df, "./data/raw/arm_df.tsv")

#Optional: cache for ProjectTemplate.
#ProjectTemplate::cache("arm_df")
