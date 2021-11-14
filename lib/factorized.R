#https://adv-r.hadley.nz/s3.html


# Factorization class -----------------------------------------------------

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
  
  return(x)
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

pca_to_factorized <- function(out, gene_names, cell_names, rank) {
  rank <- as.integer(rank)
  
  method <- "PCA"
  
  cell_mat <- out$rotation[,1:rank]
  gene_mat <- out$x[,1:rank]
  extras <- out[c("sdev", "center", "scale")]
  
  stopifnot(dim(gene_mat)[1] == length(gene_names))
  stopifnot(dim(cell_mat)[1] == length(cell_names))
  
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
fastICA_to_factorized <- function(out, gene_names, cell_names) {
  method <- "fastICA"
  
  source_mat <- out$S
  mixing_mat <- out$A
  extras <- out[names(out)[!(names(out) %in% c("A", "S"))]]
  
  rank <- ncol(source_mat)
  
  gene_mat <- source_mat
  cell_mat <- mixing_mat
  
  stopifnot(dim(gene_mat)[1] == length(gene_names))
  stopifnot(dim(cell_mat)[1] == length(cell_names))
  
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

