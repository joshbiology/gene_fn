#https://adv-r.hadley.nz/s3.html


# Factorization class -----------------------------------------------------
#Factorization class - scalar object - list represents one thing

#Constructor: efficiently creates new objects with the correct structure
new_factorized <- function(gene_mat = double(),
                           cell_mat = double(),
                           rank = integer(),
                           method = character(), 
                           gene_names = character(), 
                           cell_names = character(),
                           extras = list()) {
  
  stopifnot(is.double(gene_mat))
  stopifnot(is.double(cell_mat))
  stopifnot(is.integer(rank))
  stopifnot(is.character(method))
  stopifnot(is.character(gene_names))
  stopifnot(is.character(cell_names))
  stopifnot(is.list(extras))
  
  structure(list(gene_mat = gene_mat, 
                 cell_mat = cell_mat,
                 rank = rank), 
            class = "factorized",
            method = method,
            gene_names = gene_names,
            cell_names = cell_names,
            extras = extras)
}

# Validator ---------------------------------------------------------------


#Validator: more computationally expensive checks to ensure that the object has correct values
validate_factorized <- function(x) {
  gene_names <- attr(x, "gene_names")
  cell_names <- attr(x, "cell_names")
  method <- attr(x, "method")
  
  #Non-zero rank
  if (!(x$rank > 0)) {
    stop(
      "Factorized matrix must have non-zero rank.",
      call. = FALSE
    )
  }
  
  #Rows of both matrices should have the correct length
  if (length(gene_names) != nrow(x$gene_mat)) {
    stop(
      "The length of the gene names attribute must match the number of rows in the gene matrix.",
      call. = FALSE
    )
  }
  
  if (length(cell_names) != nrow(x$cell_mat)) {
    stop(
      "The length of the cell names attribute must match the number of rows in the cell matrix.",
      call. = FALSE
    )
  }
  
  
  #Columns should match the rank
  if (x$rank != ncol(x$gene_mat)) {
    stop(
      "The rank must match the number of columns in the gene matrix.",
      call. = FALSE
    )
  }
  
  if (x$rank != ncol(x$cell_mat)) {
    stop(
      "The rank must match the number of columns in the cell matrix.",
      call. = FALSE
    )
  }
  
  x
}

# Helper ------------------------------------------------------------------


#Helper: convenient way for others to create objects of your class
factorized <- function(gene_mat = double(),
                       cell_mat = double(),
                       rank = integer(),
                       method = character(), 
                       gene_names = character(), 
                       cell_names = character(),
                       extras = list()) {
  
  validate_factorized(new_factorized(gene_mat,
                                     cell_mat,
                                     rank,
                                     method,
                                     gene_names,
                                     cell_names,
                                     extras))
}


# Subhelpers --------------------------------------------------------------


#Subhelpers: individual conversion functions for specific factorization methods

pca_to_factorized <- function(out, rank) {
  #Right now, this is assuming an orientation in which the genes are the samples and the cell lines are the columns. 
  #If this changes, I'll have to change this function.
  #Genes IID, Cell lines as features
  rank <- as.integer(rank)
  
  method <- "PCA"
  
  cell_mat <- out$rotation[,1:rank]
  gene_mat <- out$x[,1:rank]
  extras <- out[c("sdev", "center", "scale")]
  
  gene_names <- rownames(gene_mat)
  cell_names <- rownames(cell_mat)
  
  #Clear metadata from matrices
  attr(gene_mat, "dimnames") <- NULL
  attr(cell_mat, "dimnames") <- NULL
  
  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}

#For Icasso stabilized ICA
fastICA_to_factorized <- function(out) {
  method <- "fastICA"
  source_mat <- out$S
  mixing_mat <- out$A
  extras <- out[names(out)[!(names(out) %in% c("A", "S"))]]
  
  rank <- ncol(source_mat)
  
  #Right now, all assumptions are that # cell lines < # genes. This may change in future.
  #Genes IID, Samples as features
  
  if (nrow(source_mat) > nrow(mixing_mat)) {
    gene_mat <- source_mat
    cell_mat <- mixing_mat
  }
  
  #Samples IID, Genes as features
  else {
    cell_mat <- source_mat
    gene_mat <- mixing_mat
  }
  
  gene_names <- rownames(gene_mat)
  
  cell_names <- rownames(cell_mat)
  
  #Clear metadata from matrices
  attr(gene_mat, "dimnames") <- NULL
  attr(cell_mat, "dimnames") <- NULL

  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}

#For Unstable ICA

reg_fastICA_to_factorized <- function(out) {
  method <- "fastICA_regular"
  source_mat <- out$S
  mixing_mat <- out$A %>% t()
  extras <- out[names(out)[!(names(out) %in% c("A", "S"))]]
  
  rank <- ncol(source_mat)
  
  #Right now, all assumptions are that # cell lines < # genes. This may change in future.
  #Genes IID, Samples as features
  
  if (nrow(source_mat) > nrow(mixing_mat)) {
    gene_mat <- source_mat
    cell_mat <- mixing_mat
  }
  
  #Samples IID, Genes as features
  else {
    cell_mat <- source_mat
    gene_mat <- mixing_mat
  }
  
  gene_names <- rownames(gene_mat)
  
  cell_names <- colnames(out$X)
  
  #Clear metadata from matrices
  attr(gene_mat, "dimnames") <- NULL
  attr(cell_mat, "dimnames") <- NULL
  
  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}

rica_to_factorized <- function(out, gene_names, cell_names) {
  method <- "RICA"
  rank <- out$k %>% as.integer()
  gene_mat <- out$z %>% t()
  cell_mat <- out$W
  extras <- out[c("lambda", "recon", "mse")]
  
  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}

sparsepca_to_factorized <- function(out, names, orientation = c("gene", "cell")) {
  
  stopifnot(orientation %in% c("gene", "cell"))
  
  method <- "sparsePCA"
  
  if (orientation == "cell") {
    
    gene_mat <- out$loadings
    cell_mat <- out$scores
    rank <- ncol(gene_mat)
    
    cell_names <- rownames(cell_mat)
    gene_names <- names
    #Clear metadata from matrices
    attr(cell_mat, "dimnames") <- NULL
  }
  
  else if (orientation == "gene") {
    
    gene_mat <- out$scores
    cell_mat <- out$loadings
    rank <- ncol(gene_mat)
    
    gene_names <- rownames(gene_mat)
    cell_names <- names
    #Clear metadata from matrices
    attr(gene_mat, "dimnames") <- NULL
  }
  
  extras <- out[names(out)[!(names(out) %in% c("A", "S"))]]

  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}

graphdl_to_factorized <- function(out, gene_names, cell_names) {
  method <- "graphDL"
  
  cell_mat <- out$D
  gene_mat <- t(out$X)
  extras <- out$params
  
  stopifnot(dim(gene_mat)[1] == length(gene_names))
  stopifnot(dim(cell_mat)[1] == length(cell_names))
  
  rank = as.integer(dim(gene_mat)[2])
  
  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}


ksvd_to_factorized <- function(out, gene_names, cell_names) {
  method <- "k-svd"
  
  cell_mat <- out$D
  gene_mat <- t(out$X)
  extras <- c(out$params, list(out$err))
  
  stopifnot(dim(gene_mat)[1] == length(gene_names))
  stopifnot(dim(cell_mat)[1] == length(cell_names))
  
  rank = as.integer(dim(gene_mat)[2])
  
  factorized(gene_mat,
             cell_mat,
             rank,
             method,
             gene_names,
             cell_names,
             extras)
}


# Generics ----------------------------------------------------------------


print.factorized <- function(x) {
  x <- unclass(x)
  method <- attr(x, "method")
  cat("Factorized matrix using", method, "and rank", x$rank, "\n")
  cat("# of genes:", nrow(x$gene_mat), "\n")
  cat("# of cells:", nrow(x$cell_mat), "\n")
}

get_gene_mat <- function(x) {
  out <- x$gene_mat
  rownames(out) <- attr(x, 'gene_names')
  colnames(out) <- paste("V", 1:x$rank, sep = "")
  return(out)
}

#This version scales coefficients into units of standard deviation.
get_gene_mat_stdev <- function(x) {
  out <- x$gene_mat
  rownames(out) <- attr(x, 'gene_names')
  colnames(out) <- paste("V", 1:x$rank, sep = "")
  
  out <- out/sqrt(dim(x$cell_mat)[1])
  return(out)
}

get_gene_loadings <- function(x, gene_name){
  get_gene_mat_stdev(x)[gene_name,] %>% enframe("Function", "Loading") %>% arrange(-abs(Loading))
}

get_function_loadings <- function(x, fn_num){
  get_gene_mat_stdev(x)[,fn_num] %>% enframe("Gene", "Loading") %>% arrange(-abs(Loading))
}

get_cell_mat <- function(x) {
  out <- x$cell_mat
  rownames(out) <- attr(x, 'cell_names')
  colnames(out) <- paste("V", 1:x$rank, sep = "")
  return(out)
}

#might need something here to give objects a name based on some of the parameters that it was run on.