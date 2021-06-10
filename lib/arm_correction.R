arm_correct <- function(dataset) {
  
  #Calculate a median arm score
  #Correction vector is gene median - arm median.
  
  source("./lib/gene_id_conversion.R")
  load('./cache/arm_df.RData')
  
  #If there are two locations, filter to the first one in the list.
  #Only affects a few poorly characterized and potentially artifactual orfs.
  
  arm_df <- arm_df %>% 
    unite(chr_arm, chromosome_name, arm) %>% 
    dplyr::select(entrezgene, chr_arm) %>% 
    unique() %>% 
    group_by(entrezgene) %>% 
    slice(1) %>% 
    ungroup() 
  
  
  arm_key_value <- arm_df$chr_arm
  names(arm_key_value) <- arm_df$entrezgene
  unique_arms <- unique(arm_key_value)
  
  #Arm medians by cell line
  arm_med_by_cell_line <- purrr::map(unique_arms, function(arm) {
    matrixStats::rowMedians(dataset[,arm_key_value[colnames(dataset) %>% convert_cds_to_entrez()] == arm], na.rm = T)}) 
  
  names(arm_med_by_cell_line) <- unique_arms
  
  arm_med_by_cell_line <- as.data.frame(arm_med_by_cell_line) %>% 
    as.matrix %>% 
    set_rownames(rownames(dataset)) %>% 
    set_colnames(unique_arms)
  
  
  #Global arm medians
  arm_med_global <- matrixStats::colMedians(arm_med_by_cell_line, na.rm = T) %>% 
    set_names(unique_arms)
  
  
  #Correction matrix on arm level
  arm_correction <- sweep(arm_med_by_cell_line, 2, arm_med_global) %>% 
    cbind("none" = rep(0,dim(arm_med_by_cell_line)[1]))
  
  
  #Correction matrix on gene level
  final_arms <- arm_key_value[colnames(dataset) %>% convert_cds_to_entrez()] %>% replace_na("none")

  gene_correction <- (dataset - arm_correction[,final_arms]) %>% 
    set_colnames(colnames(dataset))
  

  
  return(gene_correction)
}
