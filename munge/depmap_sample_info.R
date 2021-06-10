#sample info
#Eternal dataset
sample_info <- read_tsv(file.path(".", "data", "raw", "depmap", "eternal-depmap-datasets_v96-sample-info.tsv"))

#convert ccle_name to depmap_id

convert_ccle_to_depmap_id <- function(cell_line_names) {
  key <- pull(sample_info, 'CCLE Name')
  value <- pull(sample_info, 'DepMap_ID')
  dict <- value %>% magrittr::set_names(key)
  return(as.character(dict[cell_line_names]))
}

convert_depmap_id_to_ccle <- function(cell_line_names) {
  key <- pull(sample_info, 'DepMap_ID')
  value <- pull(sample_info, 'stripped_cell_line_name')
  dict <- value %>% magrittr::set_names(key)
  return(as.character(dict[cell_line_names]))
}
