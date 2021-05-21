#' Re-scales the category scale values and Phi after convergence of nominal model
#'
#' The function only applies to nominal models.  During estimation, a scaling
#' constraint is put on conditional variances such that they are equal to 1.
#' This function provide an alternative identification constraints that can be
#' used after the algorithm has converged. This function allow a user to tease
#' apart the strength and structure of associations.  The common alternative is
#' to constrain the sum of category scale values equals 0 and the sum of squares
#' equals 1. The phi parameters are adjusted accordingly and is now a
#' conditional covariance matrix. One item per trait needs to be selected for ID
#' constraint and this is indicated by the object anchor. This function is
#' designed for nominal models.
#'
#' @param item.log	    Log history of parameter estimates at each iteration
#' @param Phi.mat		    Final estimate of the phi matrix
#' @param anchor        Indicator of item(s) to place scaling constraint on.
#' @param item.by.trait Index for which trait item loads on
#' @param nitems        Number of items
#' @param ncat          Number of categories per item
#' @param nless         Number of unique lambdas &  unique nus
#' @param ntraits       Number of latent traits
#' @param ItemNames     Item names
#'
#' @return sNu         	Re-scaled category scale values
#' @return sPhi.mat	  	Re-scale phi matrix (conditional covariance matrix)
#'
#' @export

ScaleItem <- function(item.log, Phi.mat, anchor, item.by.trait, nitems,
                      ncat, nless, ntraits, ItemNames) {
  # --- get up-dated nus from last lines of item.log
  j1 <- 3 + nless
  j2 <- j1 + nless - 1
  AllNu <- matrix(0,nrow=nitems,ncol=nless)
  for (i in 1:nitems) {
    hist<- item.log[[i]]
    AllNu[i,1:nless] <- hist[nrow(hist),j1:j2]
  }
  Nu1 <-  -rowSums(AllNu)
  Nu1 <-  matrix(Nu1,nrow=nitems,ncol=1)
  AllNu <- cbind(Nu1,AllNu)

  # --- Get scaling constant and scale nus
  sc <- matrix(0,nrow=1,ncol=ntraits)
  for (item in 1:nitems) {
    if (anchor[item]==1) {
      sc[1,item.by.trait[,item]] <- sqrt(sum(AllNu[item,]**2))
    }
  }

  sNu.p <- matrix(NA,nrow=nitems,ncol=ncat)
  for (p in 1:ntraits) {
    for (item in 1:nitems) {
      if (item.by.trait[,item] == p) {
        sNu.p[item, ] <- AllNu[item,]/sc[1,p]
      }
    }
  }
  rownames(sNu.p) <- ItemNames

  # --- Apply scaling to Phis  (ugly but works)
  for (item in 1:nitems) {
    if (anchor[item]==1) {
      for (p in 1:ntraits) {
        if (item.by.trait[,item] == p) {
          Phi.mat[p,p] <- Phi.mat[p,p] * sc[1,p]**2
          for (q in 1:ntraits) {
            if (p != q) {
              Phi.mat[p,q] <- Phi.mat[p,q] * sc[1,p]
              Phi.mat[q,p] <- Phi.mat[p,q]
            }
          }
        }
      }
    }
  }

  results <- list(sNu=sNu.p, sPhi.mat=Phi.mat)
  return(results)
}
