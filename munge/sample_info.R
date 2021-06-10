#munge_cell_line meta

#Eternal dataset
sample_info <- load.from.taiga(data.name='depmap-a0ab', data.version=12, data.file='sample_info')

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
