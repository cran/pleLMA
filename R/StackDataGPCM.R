#' Prepares data for up-dating association parameters of LMA (GPCM) model
#'
#' Prepares data frame for input to 'mnlogit' for the stacked regression 
#' to obtain association parameters of the multidimensional LMA models
#' corresponding to the GPCM.  This function is called from 'fit.gpcm'. 
#' It generally would not run outside of 'fit.nominal' or 'ple.lma'.
#
#' @param 	Master		Current master data frame
#' @param  	item.log	Current history of iterations from which currenat 
#'                    a parameters are drawn
#' @param   starting.sv Fixed category scores
#' @param   npersons	Number of persons
#' @param   nitems		Number of items
#' @param   ncat		  Number of categories
#' @param   nless		  Number of unique lambdas
#' @param   Maxnphi		Number of estimated phi parameters
#' @param   pq.mat		Array needed for rest-total scores
#' @param   LambdaNames Names of lambdas in Master/Stacked data set
#' @param   PhiNames	Names of phi parameters
#'
#' @return  stack.data	Formats data for input to mnlogit to up-date phi parameters
#'
#' @export

StackDataGPCM <- function(Master, item.log, starting.sv, npersons, nitems,
                          ncat, nless,  Maxnphi, pq.mat, LambdaNames,
                          PhiNames) {

  NewNu <- matrix(NA,nrow=nitems,ncol=ncat)
  for (item in 1:nitems) {
    tmp <- as.matrix(item.log[[item]])
    a <- matrix(tmp[nrow(tmp), ncol(tmp)], nrow=1, ncol=1)
    NewNu[item,] <- as.vector(a) * as.matrix(starting.sv[item,])
  }
  NewNu.block <- matrix(t(NewNu), nrow=nitems*ncat, ncol=1)    # New's for each cat & item

  # --- Repeat each block for each person,
  NuItemCat <- do.call(rbind, replicate(npersons, NewNu.block, simplify=FALSE))

  # --- Needed for matrices to conform
  NuItemCat.wide <- do.call(cbind, replicate(Maxnphi,NuItemCat,simplify=FALSE))

  # --- Need nperson x nu matrix of dim (npersons x nitems).
  #  Take 1st row from every nu block
  last.lambda <- 5+nitems*nless
  PersonByNu = as.matrix(Master[seq(1, nrow(Master), nitems*ncat), (last.lambda+1):(last.lambda+nitems)])

  # Compute values to estimate the phis
  phi.col <- matrix(0, nrow=nitems*npersons, ncol=Maxnphi)
  colnames(phi.col) <- PhiNames
  irow <- 1
  for (person in 1:npersons) {
    for (item in 1:nitems) {
      for (trait.combo in 1:Maxnphi) {
        phi.col[irow,trait.combo] <- as.vector(PersonByNu[person,] %*% pq.mat[ , item , trait.combo] )
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
