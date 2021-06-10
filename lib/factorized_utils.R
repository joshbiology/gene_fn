#Utils for factorized objects

# Calculations ------------------------------------------------------------


#sparseness -- simple function calculating sparsity across a row or col of matrix.
sparseness <- function(x) {
  n <-  length(x)
  l_1 <- sum(abs(x))
  l_2 <- sqrt(sum(x^2))
  
  return((sqrt(n) - l_1/l_2)/(sqrt(n)-1))
}


#frac_active -- simple function calculating frac of activations over 1 stdev across a row or col of matrix.
frac_active <- function(x, threshold = 1) {
  #Standardize
  x <- (x - mean(x))/sd(x)
  sum(abs(x) >= threshold)/length(x)
}


#cor_spearman -- simple wrapper for spearman
cor_spearman <- function(...) {
  cor(..., method = "spearman", use = "pairwise")
}

#cosine_sim --
#https://stats.stackexchange.com/questions/31565/is-there-an-r-function-that-will-compute-the-cosine-dissimilarity-matrix
#Takes a single mat and returns cosine sim between rows
cosine_sim <- function(mat1, mat2 = mat1) {
  mat=Matrix::tcrossprod(mat1, mat2)
  t1=sqrt(apply(mat1, 1, Matrix::crossprod))
  t2=sqrt(apply(mat2, 1, Matrix::crossprod))
  mat / outer(t1,t2)
}

#cosine_dist -- 
#simple wrapper for converting cosine to dist
cosine_dist <- function(mat) {
  dst <- 1 - cosine_sim(mat)
  diag(dst) <- 0
  dst
}

#euc_dist -- 
#simple wrapper for euc. dist
euc_dist <- function(mat) {
  dist(mat) %>% as.matrix()
}

#calculate_coranking -- 
#Wrapper around co-ranking metrics using any distance function, def. cosine
calculate_coranking <- function(dr_mat, orig_mat, 
                                dr_dist_fn = cosine_dist,
                                orig_dist_fn = cosine_dist) {
  Q <- coRanking::coranking(orig_dist_fn(orig_mat),
                            dr_dist_fn(dr_mat), "dist")
  
  #Gather outputs
  lcmc <- coRanking::LCMC(Q)
  Kmax <- which.max(lcmc) %>% as.numeric()
  Q_Local <- mean(lcmc[1:Kmax])
  Q_Global <- mean(lcmc[(Kmax + 1):nrow(Q)])
  
  
  #Pass outputs
  list(LCMC = lcmc,
       K_Max = Kmax,
       Q_Local = Q_Local,
       Q_Global = Q_Global)
  
}


# Tools for Factorized objects ---------------------------------------
#Any methods should output raw data for plotting

#fgsea_on_factorized
#Input is factorized object, genes, and communities.
#I shouldn't be ingesting data frame objects, usually because there is
#cheaper function logic I can use.

fgsea_on_factorized <- function(fctr_obj, genes, communities) {
  stopifnot(is.character(genes))
  stopifnot(is.character(communities))
  
  #Enforce lengths are equal
  if (length(genes) != length(communities)) {
    stop(
      "Length of genes and communities must be the same.",
      call. = FALSE
    )
  }
  
  #Assemble input data (using factored_genes instead of gene_names to prevent confusion)
  factored_genes <- attr(fctr_obj, "gene_names")
  gene_mat <- fctr_obj$gene_mat %>% 
    set_rownames(factored_genes)
  
  #Create geneset list object, filtered for relevant genes
  index <- which(genes %in% factored_genes)
  pathways <- split(genes[index], communities[index])
  
  #Run enrichment tests
  adply(gene_mat, 2, function(data) {
    fgsea::fgsea(pathways, stats = data, minSize = 3, maxSize = 30, nperm = 1000)
  }, .id = "Factor")
}


#sparseness_on_factorized -- 
#wrapper to run on fctr_obj
sparseness_on_factorized <- function(fctr_obj) {
  #This function maps sparsity across each gene and cell line observation.
  list(Gene = adply(fctr_obj$gene_mat, 1, sparseness),
       Cell = adply(fctr_obj$cell_mat, 1, sparseness)) %>% 
    enframe("Type") %>% 
    unnest("value") %>% 
    rename(Sparseness = V1)
}


#frac_active_on_factorized -- 
#wrapper to run on fctr_obj
frac_active_on_factorized <- function(fctr_obj) {
  #This function maps the fraction of active units in a single factor.
  list(Gene = adply(fctr_obj$gene_mat, 2, frac_active, .id = "Factor"),
       Cell = adply(fctr_obj$cell_mat, 2, frac_active, .id = "Factor")) %>% 
    enframe("Type") %>% 
    unnest("value") %>% 
    rename(Frac_Active = V1)
}

#tidy_factorized -- 
#returns a tidy df from fctr_obj
tidy_factorized <- function(fctr_obj, scale = T){
  #This function maps the fraction of active units in a single factor.
  
  gene_vals <- fctr_obj$gene_mat %>% as.vector() 
  cell_vals <- fctr_obj$cell_mat %>% as.vector() 
  
  if (scale) {
    gene_vals <- gene_vals %>% scale()
    cell_vals <- cell_vals %>% scale()
  }
  
  rbind(tibble(Type = "Gene", Value = gene_vals), tibble(Type = "Cell", Value = cell_vals))
}

#recon -- 
#returns a reconstructed matrix from fctr_obj
recon <- function(fctr_obj) {
  gene_names <- attr(fctr_obj, "gene_names")
  cell_names <- attr(fctr_obj, "cell_names")
    
  Rfast::Tcrossprod(fctr_obj$cell_mat, fctr_obj$gene_mat) %>% 
    set_rownames(cell_names) %>% 
    set_colnames(gene_names)
}

#evaluate_recon_factorized -- 
#returns a df with reconstruction error calculated across rows and cols of reconstructed mat from fctr_obj
#Takes in metric functions that follow the Metrics api.
evaluate_recon_factorized <- function(fctr_obj, orig_mat, metric_fn = Metrics::mse, .parallel = F) {
  recon_mat <- recon(fctr_obj)
  
  if (dim(recon_mat)[1] != dim(orig_mat)[1] | dim(recon_mat)[2] != dim(orig_mat)[2]) {
    stop("Dimensions do not match. Rows as genes, cells as columns.",
         call. = FALSE)
  }
  
  if (.parallel) {
    gene_error <- furrr::future_map_dbl(1:nrow(recon_mat), ~metric_fn(recon_mat[.,], orig_mat[.,]))
    cell_error <- furrr::future_map_dbl(1:ncol(recon_mat), ~metric_fn(recon_mat[,.], orig_mat[,.]))
    
  }
  else {
    gene_error <- purrr::map_dbl(1:nrow(recon_mat), ~metric_fn(recon_mat[.,], orig_mat[.,]))
    cell_error <- purrr::map_dbl(1:ncol(recon_mat), ~metric_fn(recon_mat[,.], orig_mat[,.]))
  }
  
  rbind(tibble(Type = "Gene", Name = rownames(recon_mat), Metric = gene_error), 
        tibble(Type = "Cell", Name = colnames(recon_mat), Metric = cell_error))
}


#coranking_on_factorized --
#returns a coranking list with different attributes that can be used for plotting / ranking
#in current implementation, calculates gene coranking and cell line co-ranking at the same time
coranking_on_factorized <- function(fctr_obj, orig_mat) {
  #out
  list(Gene = calculate_coranking(fctr_obj$gene_mat, orig_mat),
       Cell = calculate_coranking(fctr_obj$cell_mat, t(orig_mat)))
}


#cross_cor_on_factorized
#evaluates how correlated each factor is to other factors
cross_cor_on_factorized <- function(fctr_obj) {
  
  nms <- 1:fctr_obj$rank
  
  #Prepare cor matrices
  gene_fctr_cor <- cor_spearman(fctr_obj$gene_mat)
  cell_fctr_cor <- cor_spearman(fctr_obj$cell_mat)
  
  #Eliminate redundancy
  gene_fctr_cor[lower.tri(gene_fctr_cor, diag=T)] <- NA
  cell_fctr_cor[lower.tri(cell_fctr_cor, diag=T)] <- NA
  
  
  #Set dimnames
  rownames(gene_fctr_cor) <- nms
  colnames(gene_fctr_cor) <- nms
  rownames(cell_fctr_cor) <- nms
  colnames(cell_fctr_cor) <- nms
  
  
  #Out
  list(Gene = mat.to.df(gene_fctr_cor, dat.name = "Cor"),
       Cell = mat.to.df(cell_fctr_cor, dat.name = "Cor")) %>% 
    enframe("Type") %>% 
    unnest("value")  %>% 
    na.omit()
}











