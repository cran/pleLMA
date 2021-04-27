#' Prepares data for up-dating association parameters of a multidimensional nominal LMA
#'
#' Prepares data frame for input to 'mnlogit' for the stacked regression 
#' to obtain association parameters of the multidimensional LMA models
#' corresponding to the Nominal item response model. This function is
#' called from 'fit.nominal'. It generally would not run outside of 
#' either 'fit.nominal' or 'ple.lma'.
#'
#' @param Master	  	Master data frame from which stacked data is created
#' @param	item.log	  Last row contains current scale values (item.history)
#' @param	phi.log	  	Last row contains current estimates of phi
#' @param pq.mat      Summing array to get rest scores and totals
#' @param	npersons    Number of persons
#' @param	nitems      Number of items
#' @param ncat        Number of categories per item
#' @param nless       Number of categories less 1 (i.e., unique lambdas and unique nus)
#' @param	ntraits     Number of latent traits
#' @param	Maxnphi     Number of phis to be estimated
#' @param PhiNames    Names of the Phi parameters
#' @param LambdaNames Names of lambdas that correspond to those in Master
#'
#' @return	Phi.mat		Up-dated matrix of phi parameters
#' @return	phi.log 	History of iterations for LogLike, Lambda and phi parameters

#' @export

StackData <- function(Master, item.log, phi.log, pq.mat, npersons, nitems, ncat,
                      nless, ntraits, Maxnphi, PhiNames, LambdaNames) {

  last.lambda <- 5+nitems*nless

  NewNu <- matrix(0,nrow=nitems,ncol=nless)
  for (item in 1:nitems) {
    tmp <- as.matrix(item.log[[item]])
    NewNu[item,] <- matrix(tmp[nrow(tmp), (3+nless):(2+2*nless)], nrow=1, ncol=nless)
  }
  Nu1 <- -rowSums(NewNu)
  NewNu <- cbind(Nu1, NewNu)                                  # Combine Nu1 with other nus

  NewNu.block <- matrix(t(NewNu), nrow=nitems*ncat, ncol=1)    # New's for each cat & item

  # --- Repeat each block for each person,
  NuItemCat <- do.call(rbind, replicate(npersons, NewNu.block, simplify=FALSE))

  # --- Needed for matrices to conform
  NuItemCat.wide <- do.call(cbind, replicate(Maxnphi,NuItemCat,simplify=FALSE))

  # --- Need nperson x nu matrix of dim (npersons x nitems).
  #  Take 1st row from every nu block
  PersonByNu = as.matrix(Master[seq(1, nrow(Master), nitems*ncat), (last.lambda+1):(last.lambda+nitems)])

  # Compute values to estimate the phis
  phi.col <- matrix(0, nrow=nitems*npersons, ncol=Maxnphi)
  colnames(phi.col) <- PhiNames
  irow <- 1
  for (person in 1:npersons) {
    for (item in 1:nitems) {
      for (trait.combo in 1:Maxnphi) {
        phi.col[irow, trait.combo] <- as.vector(PersonByNu[person,] %*% pq.mat[ , item, trait.combo] )
      }
      irow <- irow + 1
    }
  }

  # --- Bump to Nstack
  Phi.col <- phi.col[rep(seq_len(nrow(phi.col)), each = ncat), ]

  # --- Weight rest.scores and totals by the nus for that item and category
  Phi.col <- as.data.frame(NuItemCat.wide * Phi.col)

  stack.data <- cbind(Master[,2:5], Master[,6:last.lambda], Phi.col, Master$alt, Master$choice)
  names(stack.data) <- c("CaseID", "Item", "Category", "y", LambdaNames, PhiNames, "alt", "choice")

  return(stack.data)
}
