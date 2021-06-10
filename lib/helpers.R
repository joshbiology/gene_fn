
export_for_matlab <- function(mat, out_file) {
  #Strip row and col names
  #Correct export params only avail for write.table
  write.table(mat %>% 
                set_rownames(NULL) %>% 
                set_colnames(NULL),
              row.names = F,
              col.names = F,
              sep = ",",
              file = out_file)
  
}

normalize_cols <- function(x) {
  x %*% diag(1/sqrt(colSums(x * x)))
}


umap_to_df <- function(umap_obj, rowname) {
  umap_obj$layout %>% 
    as.data.frame() %>% 
    rownames_to_column(rowname)
  
}


edgeweight_symmetric_rank <- function(mat, rank_method = "maximal") {
  edge_rank <- mat
  diag(edge_rank) <- NA
  edge_rank <- matrixStats::rowRanks(-edge_rank)  %>%
    set_rownames(rownames(edge_rank)) %>%
    set_colnames(colnames(edge_rank))
  
  if (rank_method=="maximal") {
    edge_rank <- pmin(edge_rank, t(edge_rank), na.rm=T)
  }
  if (rank_method=="minimal") {
    edge_rank <- pmax(edge_rank, t(edge_rank), na.rm=T)
  }
  if (rank_method=="average") {
    edge_rank <- (edge_rank + t(edge_rank)) / 2
  }
  return(edge_rank)
}

# Dual graph regularized dictionary learning (DGRDL) utils ------------------------
import_graphdl <- function(matlab_obj) {
  
  index <- names(matlab_obj)
  params <- matlab_obj[!(index %in% c('D', 'X'))]
  D <- matlab_obj$D
  X <- matlab_obj$X %>% as.matrix()
  
  
  list(params = params,
       D = D,
       X = X)
}

