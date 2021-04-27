#' Re-scales the category scale values and Phi after convergence of the nominal model
#'
#' This auxillary function only applies to nominal models after estimating the
#' parameters of the model.  During estimation that scaling identification
#' constraint is put on conditional variances (i.e., phi_mm) such that they
#' equal 1. This function provide an alternative identification constraint
#' after the algorithm has converged. This function allow a user to tease
#' apart the strength and structure of associations. The alternative
#' scaling identification constraint is to set the sum of category scale
#' values equals 0 and the sum of squares equal to 1. The phi parameters
#' are adjusted accordingly. Only one item per trait should be selected for
#' the identification constraint and this is indicated by the object anchor.
#'
#' @param model.fit     Model object for a nominal model
#' @param anchor        Indicator of item(s) to place scaling constraint on.
#'
#' @return sNu         	Re-scaled category scale values
#' @return sPhi.mat	  	Re-scale phi matrix (conditional covariance matrix)
#'
#' @examples
#' 
#' \donttest{
#' #--- 3 items from depression, anxiety, and stress scales
#' #    for 250 cases out of possible 1000
#' data(dass)
#' inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#' inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#' inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#' #--- nominal response model
#' n1 <- ple.lma(inData, model.type="nominal",inItemTraitAdj,inTraitAdj,tol=1e-03)
#' anchor <- c(1,0,0,0,0,0,0,0,0)
#' reScaleItem(model.fit=n1, anchor)
#'
#' #--- Multidimensional models
#' inTraitAdj  <- matrix(1, nrow=3, ncol=3)
#' dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
#' anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
#' stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
#' das <- list(dpress, anxiety, stress)
#' inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])
#'
#' n3 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
#' anchor <- c(1,0,0, 0,1,0, 1,0,0)
#' reScaleItem(model.fit=n3, anchor)
#'
#' }
#'
#' @export

reScaleItem <- function(model.fit, anchor) {

   item.log      <- model.fit$item.log
   Phi.mat       <- model.fit$Phi.mat
   item.by.trait <- model.fit$item.by.trait
   nitems        <- model.fit$nitems
   ncat          <- model.fit$ncat
   nless         <- model.fit$nless
   ntraits       <- model.fit$ntraits
   ItemNames     <- model.fit$ItemNames

  # check anchor
  check <- matrix( 0, nrow=1, ncol=ntraits)
  for (i in 1:nitems){
    for (k in 1:ntraits){
      if (item.by.trait[i]==k) { check[1,k] <- check[1,k] + anchor[i]
      }
    }
  }
  check0 <- rep(1, ntraits)

  if (sum(check0-check) != 0 ) {
    stop("There is an error in anchor. One item per trait should be specified.")
  }


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
