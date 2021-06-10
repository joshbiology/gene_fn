#adapted from MineICA::clusterFastICARuns
#bypasses need for icasso stablilization code in Matlab
#https://rdrr.io/bioc/MineICA/src/R/runAn.R
icasso_ica <- function (X, nbComp, nbIt = 100, alg.type = c("deflation", "parallel"), 
                        fun = c("logcosh", "exp"), maxit = 500, tol = 10^-6, funClus = c("hclust", 
                                                                                         "agnes", "pam", "kmeans"), row.norm = FALSE, bootstrap = FALSE, 
                        ...) 
{
  alg.type <- match.arg(alg.type)
  funClus <- match.arg(funClus)
  
  require(doParallel)
  message(paste("FastICA iteration 1"))
  resICA <- fastICA::fastICA(X = X, n.comp = nbComp, alg.typ = alg.type, 
                             fun = fun, maxit = maxit, tol = tol, row.norm = row.norm)
  whit <- resICA$K
  print("solve 1")
  
  dewhit <- solve(t(resICA$K) %*% resICA$K) %*% t(resICA$K)
  it <- clus <- NULL
  allW <- foreach(it = 2:nbIt, .combine = cbind) %dopar% {
    message(paste("FastICA iteration "), it)
    if (bootstrap) 
      Xbis <- X[sample(1:ncol(X), replace = TRUE), ]
    else Xbis <- X
    res <- fastICA::fastICA(X = Xbis, n.comp = nbComp, alg.typ = alg.type, 
                   fun = fun, maxit = maxit, tol = tol, row.norm = row.norm)
    res$W
  }
  allW <- cbind(resICA$W, allW)
  allWdewith <- t(allW) %*% dewhit
  allWdewith <- (apply(allWdewith, 2, scale))
  sim <- abs(cor(t(allWdewith)))
  dsim <- 1 - sim
  centrotypes <- c()
  switch(funClus, hclust = {
    resClus <- hclust(d = as.dist(dsim), ...)
    partition <- cutree(resClus, k = nbComp)
  }, agnes = {
    resClus <- agnes(x = dsim, diss = TRUE, ...)
    partition <- cutree(as.hclust(resClus), k = nbComp)
  }, pam = {
    resClus <- pam(x = dsim, diss = TRUE, k = nbComp, keep.diss = FALSE, 
                   ...)
    partition <- resClus$clustering
    centrotypes <- resClus$medoids
  }, kmeans = {
    resClus <- kmeans(x = allW, centers = nbComp, ...)
    partition <- resClus$cluster
  })
  Iq <- foreach(clus = unique(partition), .combine = c) %do% 
    {
      indC <- which(partition == clus)
      if (length(indC) > 1) {
        if (funClus != "pam") 
          centrotypes <- c(centrotypes, indC[which.max(apply(sim[indC, 
                                                                 indC], 1, sum))])
        internalSim <- mean(sim[indC, indC])
      }
      else {
        if (funClus != "pam") 
          centrotypes <- c(centrotypes, indC)
        internalSim <- sim[indC, indC]
      }
      externalSim <- mean(sim[indC, setdiff(1:ncol(sim), 
                                            indC)])
      iq <- internalSim - externalSim
      return(iq)
    }
  W <- whit %*% allW[, centrotypes]
  print("solve 2")
  A <- solve(t(W) %*% W) %*% t(W)
  S <- X %*% W
  rownames(S) <- rownames(X)
  colnames(A) <- colnames(X)
  orderIq <- order(Iq, decreasing = TRUE)
  return(list(A = t(A)[, orderIq], S = S[, orderIq], W = W[, 
                                                           orderIq], Iq = Iq[orderIq]))
}



