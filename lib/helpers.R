#create_output_folder -- 
#simple function that checks for the existence of a folder before creating it.
create_output_folder <- function(out_path) {
  ifelse(!dir.exists(out_path), dir.create(out_path, recursive = TRUE), FALSE)
}

#export_for_matlab -- 
#strips attributes from matrix before saving it at out_file.
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

#normalize_cols -- 
normalize_cols <- function(x) {
  x %*% diag(1/sqrt(colSums(x * x)))
}

#umap_to_df -- 
#Simple conversion tool for umap objects returned by the umap package.
umap_to_df <- function(umap_obj, rowname) {
  umap_obj$layout %>% 
    as.data.frame() %>% 
    rownames_to_column(rowname)
  
}

#edgeweight_symmetric_rank -- 
#Takes a matrix as input and outputs a matrix of same dimensions, 
#but with values that are replaced by symmetric ranked edgeweights.
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


#generate_graph -- 
#Takes a matrix as input and outputs a graph with edges between nearest neighbors (cosine similarity).

generate_graph <- function(mat, rank = 5) {
  
  tmp <- cosine_sim(mat) %>%
    edgeweight_symmetric_rank
  
  tmp[tmp >rank] <- NA
  
  tmp[upper.tri(tmp)] <- NA
  
  tmp[!is.na(tmp)] <- 1
  
  tmp2 <- igraph::graph_from_adjacency_matrix(tmp, mode = "undirected", weighted = NULL)
  
  return((tmp2))
  
}

# Dual graph regularized dictionary learning (DGRDL) utils ------------------------
#import_graphdl -- converts matlab object into list.
import_graphdl <- function(matlab_obj) {
  
  index <- names(matlab_obj)
  params <- matlab_obj[!(index %in% c('D', 'X'))]
  D <- matlab_obj$D
  X <- matlab_obj$X %>% as.matrix()
  
  
  list(params = params,
       D = D,
       X = X)
}

#extract_atoms --returns the genes and columns involved.
extract_atoms <- function(mat, gene, loading_threshold = 8) {
  # set col names
  colnames(mat) <- paste("V", 1:ncol(mat), sep = "")
  
  gene_loadings <- mat[gene, ]
  col_index <- abs(gene_loadings) > 0
  
  row_tmp <- rowSums(abs(mat[,col_index]))
  row_index <- row_tmp > loading_threshold
  
  
  return(mat[row_index, col_index])
}

