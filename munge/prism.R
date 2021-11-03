# Prism -------------------------------------------------------------------
#primary-screen-e5c7
prism_lfc <- read_csv(file.path(".", "data", "raw", "prism", "primary-screen-public-tentative_v11-primary-replicate-collapsed-logfold-change.csv")) %>%
  column_to_rownames("X1") %>%
  as.matrix()
prism_meta <- read_tsv(file.path(".", "data", "raw", "prism", "primary-screen-public-tentative_v11-primary-replicate-collapsed-treatment-info.tsv"))

#secondary-screen-0854
#PR500 - adherent cell line panel
#PR300 - mixture of adherent and suspension with some overlap to PR500 for comparison
#we keep PR500 for analysis
prism_secondary_meta <- read_tsv(file.path(".", "data", "raw", "prism","secondary-screen-public-tentative_v18-secondary-replicate-collapsed-treatment-info.tsv"))


prism_secondary_lfc <- read_csv(file.path(".", "data", "raw", "prism",
                                          "secondary-screen-public-tentative_v18-secondary-replicate-collapsed-logfold-change.csv")) %>%
  column_to_rownames("X1") %>%
  as.matrix()

prism_secondary_rownames <- read_tsv(file.path(".", "data", "raw", "prism", "secondary-screen-public-tentative_v18-cell-line-info.tsv")) %>%
  dplyr::filter(str_sub(row_name, 1, 5) == "PR500", !is.na(depmap_id), row_name %in% rownames(prism_secondary_lfc)) %>% 
  group_by(depmap_id) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()


#Filter PR500
prism_secondary_lfc <- prism_secondary_lfc[prism_secondary_rownames$row_name,] %>%
  set_rownames(prism_secondary_rownames$depmap_id)

#Sanity check
length(rownames(prism_secondary_lfc))
length(unique(rownames(prism_secondary_lfc)))


#UMAP raw data --------------------------------------------------------------
#https://depmap.org/repurposing/
#User note - a compound can have many MOA's, resulting in multiple rows. Keep that in mind when subsetting.
prism_umap_annot <- read_csv(file.path(".", "data", "raw", "prism", "prism_repurposing_umap.csv"))
