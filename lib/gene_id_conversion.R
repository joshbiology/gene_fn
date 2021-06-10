#Gene ID's

#HUGO protein coding 9-18-2019 download
#Note: the current 19Q4 release is behind on the symbol names. To maintain consistency, I should always be parsing to the entrez names.
hugo <- read_tsv("./data/raw/hugo-protein-coding_gene-9-18-2019.txt") %>% 
  mutate(tmp = paste("(", entrez_id, ")", sep = "")) %>% 
  unite_("cds_id",  c("symbol", "tmp"), remove = F, sep = " ") %>% 
  dplyr::select(-tmp)

convert_genes <- function(genes, 
                          from = c("entrez_id", "hgnc_id", "symbol", "cds_id", "ensembl_gene_id", "vega_id", "ucsc_id", "refseq_accession", "uniprot_ids", "omim_id"),
                          to = c( "symbol", "hgnc_id","entrez_id", "ensembl_gene_id", "vega_id", "ucsc_id", "refseq_accession", "uniprot_ids", "omim_id",
                                 "cds_id")) {
  
  from <- match.arg(from)
  to <- match.arg(to)
  #Relevant columns: hgnc_id, symbol, entrez_id, ensembl_gene_id, vega_id, ucsc_id, refseq_accession, uniprot_ids, omim_id
  key <- pull(hugo, from)
  
  value <- pull(hugo, to)
  
  dict <- value %>% magrittr::set_names(key)
  
  return(as.character(dict[genes]))
  #Use key-value pairs to convert from different columns in the hugo database
}

convert_cds_to_entrez <- function(genes) {
  #A separate function is required as the cds_id is not stable. 
  #Expected input is a set of cds_id names.
  str_split(genes, " ", simplify = T)[,2] %>% str_sub(2, -2)
}

convert_cds_to_symbol <- function(genes) {
  #A separate function is required as the cds_id is not stable. 
  #Expected input is a set of cds_id names.
  str_split(genes, " ", simplify = T)[,1]
}

convert_genes_mygeneinfo <- function(genes, to = c("entrezgene", "symbol")) {
  to <- match.arg(to)
  
  dat <- mygene::queryMany(genes, 
                    scopes = "entrezgene,symbol,alias, ensembl.gene", 
                    fields = "entrezgene,name,symbol,taxid,type_of_gene", 
                    species="human", 
                    size = 1,
                    entrezonly = T) %>% 
    as_tibble()
  
  key <- pull(dat, query)
  
  value <- pull(dat, to)
  
  dict <- value %>% magrittr::set_names(key)
  
  return(dict[genes] %>% as.character())
}

to_entrez <- function(symbol) {
  convert_genes(symbol, "symbol", "entrez_id")
}
