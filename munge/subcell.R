#Subcellular localization
#BioID gene map: https://humancellmap.org
bioid_loc <- read_tsv("./data/raw/humancellmap/preys-latest.txt") %>%
  rename(entrezgene = Entrez, Location = `MMF localization`, NMF_Rank = `NMF rank`) %>%
  dplyr::mutate(entrezgene = as.character(entrezgene),
         CDS_ID = convert_genes(entrezgene, "entrez_id", "cds_id"))

#BioID NMF localization profile. Each row sums to 1.
bioid_nmf <- read_excel("./data/raw/humancellmap/media-8.xlsx", sheet = 2, .name_repair = "universal") %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  aaply(1, function(x)x/sum(x))

bioid_meta <- read_excel("./data/raw/humancellmap/media-9.xlsx", sheet = 2, col_names = c("rank", "GO_Terms", "Location", "GO_IDS"), skip = 1)%>%
  dplyr::mutate(NMF_Column_Name = paste("rank", rank, sep = "."))
