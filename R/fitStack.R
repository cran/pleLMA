#' Up-dates association parameters of the nominal model
#'
#' Discrete choice model (conditional multinominal logistic regression model)
#' is fit to stacked data to up-date the matrix of association parameters of
#' the LMA that corresponds to the nominal item response model. This is a
#' function internal to 'fit.nominal' and is used for multi-dimensional models.
#' The function is similar to 'fit.StackGPCM'. This function is unlikely to
#' be run outside of 'fit.nominal' or 'ple.lma'.
#'
#' @param Master	    Master data set from which stacked data is created
#' @param item.log  	Last row contains current scale values (item.history)
#' @param	phi.log		  Last row contains current estimates of phi
#' @param fstack		  Formula for stacked regression
#' @param TraitByTrait inTraitAdj matrix
#' @param pq.mat      Summing array to get rest scores and totals
#' @param	npersons    Number of persons
#' @param	nitems      Number of items
#' @param ncat        Number of categories per item
#' @param nless       Number of categories less 1 (unique lambdas & unique nus)
#' @param	ntraits     Number of latent traits
#' @param	Maxnphi     Number of phis to be estimated
#' @param PhiNames    Names of the Phi parameters
#' @param LambdaNames Names of lambdas that correspond to those in Master
#'
#' @return	Phi.mat		 Matrix of up-dated estimates of assocation (phi) parameters
#' @return	phi.log 	 History of iterations log likelihood and estimates of
#'                       lambda and phi parameters

#' @export
FitStack <- function(Master, item.log, phi.log, fstack, TraitByTrait, pq.mat,
                     npersons, nitems, ncat, nless, ntraits, Maxnphi, PhiNames,
                     LambdaNames) {

  ########--- First stack the regressions in one data set
  last.lambda <- 5+nitems*nless

  NewNu <- matrix(0,nrow=nitems,ncol=nless)
  for (item in 1:nitems) {
    tmp <- as.matrix(item.log[[item]])
    NewNu[item,] <- matrix(tmp[nrow(tmp), (3+nless):(2+2*nless)], nrow=1,
                           ncol=nless)
  }
  Nu1 <- -rowSums(NewNu)

  # --- Combine Nu1 with other nus
  NewNu <- cbind(Nu1, NewNu)

  # New's for each cat & item
  NewNu.block <- matrix(t(NewNu), nrow=nitems*ncat, ncol=1)

  # --- Repeat each block for each person,
  NuItemCat <- do.call(rbind, replicate(npersons, NewNu.block, simplify=FALSE))

  # --- Needed for matrices to conform
  NuItemCat.wide <- do.call(cbind, replicate(Maxnphi,NuItemCat,simplify=FALSE))

  # --- Need nperson x nu matrix of dim (npersons x nitems).
  #  Take 1st row from every nu block
  PersonByNu = as.matrix(Master[seq(1, nrow(Master), nitems*ncat),
                                (last.lambda+1):(last.lambda+nitems)])

  # Compute values to estimate the phis
  phi.col <- matrix(0, nrow=nitems*npersons, ncol=Maxnphi)
  colnames(phi.col) <- PhiNames
  irow <- 1
  for (person in 1:npersons) {
    for (item in 1:nitems) {
      for (trait.combo in 1:Maxnphi) {
        phi.col[irow, trait.combo] <-
          as.vector(PersonByNu[person,] %*% pq.mat[ , item, trait.combo] )
      }
      irow <- irow + 1
    }
  }

  # --- Bump to Nstack
  Phi.col <- phi.col[rep(seq_len(nrow(phi.col)), each = ncat), ]

  # --- Weight rest.scores and totals by the nus for that item and category
  Phi.col <- as.data.frame(NuItemCat.wide * Phi.col)

  stack.data <- cbind(Master[,2:5], Master[,6:last.lambda], Phi.col, Master$alt,
                      Master$choice)
  names(stack.data) <- c("CaseID", "Item", "Category", "y", LambdaNames,
                         PhiNames, "alt", "choice")

  # ---Fit stacked data
  stacked <- dfidx::dfidx(stack.data, choice="y", idx=c("CaseID","alt"))
  model.fit <- mlogit::mlogit(fstack, stacked)

  # --- Prepare output
  parms <- model.fit$coefficients[1:(nitems*nless + Maxnphi)]
  logLike <- as.numeric(model.fit$logLik)
  phi.log <- rbind(phi.log, c(logLike, parms))

   #--- matrix of phi parameters
  phi.est <- parms[(nitems*nless+1):(nitems*nless+Maxnphi)]
  if (ntraits >1) {
    Phi.mat <- matrix(0, nrow = ntraits, ncol = ntraits)
    Phi.mat[lower.tri(Phi.mat, diag = TRUE)] <- phi.est
    Phi.mat <- t(Phi.mat)+Phi.mat - diag(diag(t(Phi.mat)))
  } else {
    Phi.mat <- phi.est
  }

  results <- list(Phi.mat=Phi.mat, phi.log=phi.log)
  return(results)
}
