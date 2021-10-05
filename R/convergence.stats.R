#' Computes statistics to assess convergence of the nominal model
#'
#' For the nominal model, convergence statistics are computed for each item, as
#' well as the algorithm as a whole. The main argument is the history or log
#' from fitting item regressions. The convergence statistics are the differences
#' between current values of the log likelihoods and item parameter estimates
#' and those from the previous iteration. The maximum over item of the
#' differences of the log likelihood values is used to determine convergence
#' of the pseudo-likelihood algorithm.  This function is used internally, but
#' it can also be used after fitting a model to examine how many iterations
#' are required for parameter estimates to get close to the final values
#' and whether any item parameters are still changing.
#'
#' @param  item.log      Iteration history of items' log likelihoods and parameter estimates
#' @param  nitems		 Number of items
#' @param  nless         Number of unique marginal terms (i.e., lambdas) and
#'                        unique category scale values(i.e., nus)
#' @param  LambdaName    Names of lambdas in item regressions
#' @param  NuName        Names of nu in item regressions
#'
#' @return diff.last   	 Differences between item loglikes & parameters on last two iterations
#' @return criterion.loglike  Maximum over items of the absolute value of LogLike differences
#' @return criterion.items  Sum of item differences of item parameters
#'
#'
#' @examples
#' \donttest{
#'  # 9 items from dass data for 250 cases
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#' #---   input for uni-dimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#'  #--- Uni-dimensional Nominal Model
#'  n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj,inTraitAdj, tol=1e-02)
#'
#' # Since this function in internal to fit.nominal, need to also run
#'  s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj)
#'
#'  convergence.stats(n1$item.log, n1$nitems, n1$nless, s$LambdaName, s$NuName)
#'
#' #--- Multidimensional models
#' #--- re-define inTraitAdj and inItemTraitAdj for 3 dimensional models
#'   inTraitAdj  <- matrix(1, nrow=3, ncol=3)
#'
#'   dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
#'   anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
#'   stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
#'   das <- list(dpress, anxiety, stress)
#'   inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])
#'
#' #--- 3 dimensional nominal
#'   n3 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
#'   s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj)
#'   convergence.stats(n3$item.log, n3$nitems, n3$nless, s$LambdaName, s$NuName)
#' }
#'
#' @export
##################################################################
convergence.stats <- function(item.log, nitems, nless, LambdaName, NuName) {

	nnn <- 3 + 2*nless
	diff.last <- matrix(0, nrow=nitems, ncol=nnn)
	for (item in 1:nitems) {
		history <- as.matrix(item.log[[item]])
		tmp <- apply(history, 2, diff)
		diff.last[item,1:nnn] <- c(item, tmp[nrow(tmp),])
	}
	diff.last <- as.data.frame(diff.last)
	diff.last <- diff.last[,c(1,3:ncol(diff.last))]

  names(diff.last) <- c("Item", "LogLike", LambdaName, NuName)

# maxiumn difference of loglike over items
  criterion.loglike <- max(diff.last$LogLike)

# sum overs all item differences.
  criterions.items  <- sum(diff.last[,3:(2+2*nless)])

  results <- c(diff.last=diff.last, criterion.loglike=criterion.loglike,
             criterion.items=criterions.items )
return(results)


}

