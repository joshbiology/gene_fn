#Genesets
parse_genesets <- function(named_paths) {
  named_paths %>%
    map(read_tsv) %>%
    enframe(name = "Process") %>%
    unnest() %>%
    dplyr::mutate(entrezgene = convert_genes(Gene, from = "symbol", to = "entrez_id")) %>%
    distinct() %>%
    rename(symbol = Gene)
}

default_genesets <- list.files("./data/interim/genesets", full.names = T, pattern = ".txt") %>%
  set_names(list.files("./data/interim/genesets", full.names = F, pattern = ".txt") %>% word(1, sep = "\\."))

genesets_df <-  default_genesets %>%
  parse_genesets()

#From HUGO gene families, downloaded 10-29-2019
#https://biomart.genenames.org

hugo_df <- read_tsv("./data/raw/hugo_gene_families.txt") %>%
  select(2, 3, 8) %>%
  set_colnames(c("Geneset", "Root_Symbol", "symbol"))

#Individual features
ints_genesets <- genesets_df %>%
  dplyr::filter(Geneset %in% c("INTS10-13-14-C7",
                        "INTS9-BRAT1-WDR73",
                        "Z3",
                        "Integrator"),
         Process == "top_hits") %>%
  select(Geneset, symbol)

dynein_genesets <- genesets_df %>%
  dplyr::filter(Geneset %in% c(#"Dynein",
    "Dynein_activators",
    #"Spindly SPDL1",
    "Dynactin",
    "SKA")) %>%
  dplyr::filter(!is.na(symbol)) %>%
  select(Geneset, symbol)

pausing_genesets <- rbind(genesets_df %>%
                            dplyr::filter(Geneset %in% c("DSIF", "PAF")) %>%
                            select(Geneset, symbol),
                          hugo_df %>%
                            dplyr::filter(Geneset %in% c("Negative elongation factor complex members",
                                                  "Super elongation complex")) %>%
                            select(Geneset, symbol)) %>%
  distinct()


mediator_genesets <- genesets_df %>%
  dplyr::filter(Process == "mediator_head_middle_ckm") %>%
  select(Geneset, symbol)


baf_genesets <- genesets_df %>%
  dplyr::filter(Process == "baf") %>%
  dplyr::filter(symbol %in% c("ARID1A", "SMARCC1", "SMARCB1", "SMARCE1", "SMARCD1", "BRD9", "BICRA", "BRD7", "ARID2", "PBRM1"))





# Hallmark Genesets -------------------------------------------------------
#Iorio SLAP paper

iorio_hallmark <- readxl::read_excel('./data/raw/iorio_hallmark_pathays.xlsx') %>%
  select(-1) %>%
  rename(Pathway = 2, Source = 3) %>%
  separate_rows(Genes)

# Signalling Genesets -----------------------------------------------------
#2018 TCGA paper

pathway_names <- readxl::excel_sheets('./data/raw/tcga_2018_signalling.xlsx')
tmp1 <- readxl::read_excel('./data/raw/tcga_2018_signalling.xlsx', sheet = 1, skip = 2) %>% select(1,3)
tmp2 <- map(pathway_names[2:10], ~readxl::read_excel('./data/raw/tcga_2018_signalling.xlsx', sheet = .) %>% select(1,3))


tcga_signalling <- c(list(tmp1), tmp2) %>%
  set_names(pathway_names[1:10]) %>%
  ldply(.id = "Pathway")

rm(pathway_names)
rm(tmp1)
rm(tmp2)
