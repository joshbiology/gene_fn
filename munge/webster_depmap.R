#munge reference webster run for depmap
#requires ./munge/helpers.R to be loaded


load("./cache/avana_19q4_webster.RData")
mat_path <- file.path(".", "data", "interim", "webster_depmap_freeze_2021-05-24-depmap_deep;K_param=220;T_param=4;alpha=0.2;beta=0.6;iternum=20;seed=1;num_neighbor_gene=5;num_neighbor_cl=5.mat")
webster_depmap <- graphdl_to_factorized(import_graphdl(R.matlab::readMat(mat_path)), 
                      colnames(avana_19q4_webster), 
                      rownames(avana_19q4_webster))


# Flip axes consistently --------------------------------------------------

#Webster solutions are sign invariant. In practice, tend to flip them so that the highest loadest gene is positive.

direction <- abs(matrixStats::colMaxs(get_gene_mat(webster_depmap))) - abs(matrixStats::colMins(get_gene_mat(webster_depmap)))

#Flip axis of ATAC factor for positive loadings.
direction[196] <- -direction[196] 

direction_mat <- diag(direction/abs(direction))

webster_depmap$gene_mat <- webster_depmap$gene_mat %*% direction_mat
webster_depmap$cell_mat <- webster_depmap$cell_mat %*% direction_mat



# Summary metrics ---------------------------------------------------------

recon_depmap <- recon(webster_depmap)
gene_recon_webster <- tibble(Gene = colnames(avana_19q4_webster),
                             Recon_Pearson = map_dbl(1:ncol(avana_19q4_webster), ~cor(avana_19q4_webster[,.], recon_depmap[,.])))

