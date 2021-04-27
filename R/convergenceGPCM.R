#' Computes statistics to assess convergence for generalized partial credit models
#'
#' For the generalized partial credit model, convergence statistics are computed
#' for each item, as well as the algorithm as a whole. The convergence statistics
#' are the differences between current values of the log likelihoods and item
#' parameter estimates and those from the previous iteration. The maximum over
#' items' differences of the log likelihood values is used to determine
#' convergence of the pseudo-likelihood algorithm. This function is used
#' internally, but it can also be used after fitting a model to examine how
#' many iterations are required for parameter estimates to get close to the
#' final values and whether any item parameters are still changing.
#'
#' @param  item.log	  Iteration history of items' log likelihoods and parameter estimates
#' @param  nitems			Number of items
#' @param  ncat				Number of categories
#' @param  nless      Number of unique lambdas
#' @param  LambdaName Names of lambdas in formula used to fit item regressions (internal to fit_gpcm)
#'
#' @return diff.last   	       Differences between item's log likelihoods and parameters for each item
#' @return criterion.loglike   Maximum overs item of the absolute value of differences between item LogLike values
#' @return criterions.items    Sum of item differences of their parameters
#'
#' @examples
#' 
#' # 9 items from dass data for 250 cases
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#'  #---   input for uni-dimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'  
#' \donttest{
#'  #--- fit Unidiemsional gpcm Model
#'  g1<- ple.lma(inData, model.type="gpcm",inItemTraitAdj,inTraitAdj, tol= 1e-03)
#'
#'  # Since convergenceGPCM is internal to fit.gpcm, need to get 'Lambdaname'
#'  s <- set.up(inData, model.type='gpcm', inTraitAdj, inItemTraitAdj)
#'
#'  convergenceGPCM(g1$item.log, g1$nitems, g1$ncat, g1$nless, s$LambdaName)
#' 
#' #--- Multidimensional models
#' #--- re-define inTraitAdj and inItemTraitAdj for 3 dimensional models
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
#'  dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
#'  anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
#'  stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
#'  das <- list(dpress, anxiety, stress)
#'  inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])
#'
#' #--- 3 dimensional gpcm
#'  g3 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1e-03)
#'  s <- set.up(inData, model.type='gpcm', inTraitAdj, inItemTraitAdj)
#'  convergenceGPCM(g1$item.log, g1$nitems, g1$ncat, g1$nless, s$LambdaName)
#' }
#'
#' @export
convergenceGPCM <- function(item.log, nitems, ncat, nless, LambdaName) {
	nnn <- 2 + ncat
	diff.last <- matrix(, nrow=nitems,ncol=nnn)
	for (item in 1:nitems) {
		hist <- as.matrix(item.log[[item]])
		hist <- hist[,2:(ncat+2)]
		tmp <- apply(hist, 2, diff)
		diff.last[item,] <- c(item,tmp[nrow(tmp),])
	}
	diff.last <- as.data.frame(diff.last)

	names(diff.last) <- c("Item", "LogLike", LambdaName, "alpha")

	criterion.loglike <- max(abs(diff.last$LogLike)  )      # maximum loglike over items
	criterions.items  <- sum(abs(diff.last[,3:(2+ncat)]))   # sum overs all item differences.

	results <- c(diff.last = diff.last,
				criterion.loglike = criterion.loglike,
				criterions.items=criterions.items )
	return(results)
}
