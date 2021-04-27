#' Fits an LMA using fixed category scores
#'
#' The LMA model with fixed category scores is fit by this function and
#' the model corresponds to models in the Rasch family of item response
#' models.  The category scores can be set by either the user or the
#' package defaults. The default category scores are equally spaced, sum to
#' zero, and sum of squares equal 1. Scores can be set by user by
#' specifying them in the item by category matrix of 'starting.sv'.  The
#' pseudo-likelihood algorithm only runs a single stacked regression.
#' This functionis called from' ple.lma' but can also be run outside
#' of the main wrapper function.
#'
#' @param	 Master        Master data set in long format
#' @param  npersons      Number of persons
#' @param  nitems		     Number of items
#' @param  ncat		       Number of categories
#' @param  nless         Number of unique Lambdas (i.e., ncat-1)
#' @param  Maxnphi       Number of association parameters
#' @param  pq.mat        One dimensional array to compute rest-scores
#' @param  starting.sv   Fixed category scores
#' @param  LambdaNames   Names of lambda paramters in Master and formula
#'                         for stacked regression
#' @param  PhiNames      Names of association parameters
#' @param  ItemNames     Names of items
#' @param  LambdaName    Names of lambdas used in output
#'
#' @return estimates 	  An item by parameter matrix of the maximum of the log likelihood,
#'                        estimated item parameters (i.e., Lambdas), and the values of
#'                        the fixed category scores.
#' @return fstack	    	Formula for stacked regression
#' @return phi.mnlogit	Results from mnlogit for stacked regression
#' @return estimates    An item x parameter estimate matrix and fixed category scores used
#' @return Phi.mat      Estimated phi parameters
#' @return mlpl.phi     Value of maximum of log pseudo-likelihood function from the stacked regression
#' @return AIC          Akaike information criterion for pseudo-likelihood (smaller is better)
#' @return BIC          Bayesian information criterion for pseudo-likelihood (smaller is better)
#'
#' @examples
#'  #---  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'  #--- unidimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#'  s <- set.up(inData, model.type='rasch', inTraitAdj, inItemTraitAdj)
#'
#'  r <- fit.rasch(s$Master, s$npersons, s$nitems, s$ncat, s$nless, s$Maxnphi,
#'           s$pq.mat, s$starting.sv, s$LambdaNames, s$PhiNames, s$ItemNames,
#'           s$LambdaName)
#'
#' @export
fit.rasch <- function(Master, npersons, nitems, ncat, nless, Maxnphi, pq.mat,
                      starting.sv, LambdaNames, PhiNames, ItemNames, LambdaName) {

  # category scores
  NewNu.block <- matrix(t(starting.sv), nrow=nitems*ncat, ncol=1)

  # --- Repeat each block for each person,
  NuItemCat <- do.call(rbind, replicate(npersons, NewNu.block, simplify=FALSE))

  # --- So that matrices conform
  NuItemCat.wide <- do.call(cbind, replicate(Maxnphi,NuItemCat,simplify=FALSE))

  # --- Need nperson x nu matrix of dim (npersons x nitems).
  #  Take 1st row from every nu block
  last.lambda <- 5+nitems*nless
  PersonByNu = as.matrix(Master[seq(1, nrow(Master), nitems*ncat),
                                (last.lambda+1):(last.lambda+nitems)])

  # Compute values to estimate the phis
  phi.col <- matrix(0, nrow=nitems*npersons, ncol=Maxnphi)
  colnames(phi.col) <- PhiNames
  irow <- 1
  for (person in 1:npersons) {
    for (item in 1:nitems) {
      for (trait.combo in 1:Maxnphi) {
        phi.col[irow,trait.combo] <-
          as.vector(PersonByNu[person,] %*% pq.mat[ , item , trait.combo] )
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

  # Formulas
  # --- formula for stacked regression, coefficients are lambdas and phis
  xstack.names <- c(LambdaNames,PhiNames)
  fstack <- stats::as.formula(paste("y ~",
                                    paste(xstack.names,collapse="+"),"| 0 | 0",sep=" "))

  # Fit model
  #--- fit the model
  phi.mnlogit <- mnlogit::mnlogit(fstack, stack.data, choiceVar="Category")

  #---- Prepare output

  # --- table of lambdas
  parms <-  phi.mnlogit$coefficients
  estimates <- matrix(NA,nitems,nless)
  x <- matrix(NA,nrow=1,ncol=ncat)
  i <- 1
  for (item in 1:nitems) {
    for (cat in 1:nless) {
      estimates[item,cat] <- parms[i]
      x[cat] <- paste("x", cat, sep="")
      i <- i +1
    }
  }
  x[ncat] <- paste("x", ncat, sep="")

  estimates <- cbind(-rowSums(estimates), estimates, starting.sv)
  rownames(estimates) <- ItemNames
  colnames(estimates) <- c("lam1", LambdaName, x)

  #--- matrix of phi parameters
  phi.est <- phi.mnlogit$coefficients[(nitems*nless+1):(nitems*nless+Maxnphi)]
  Phi.mat <- phi.est

  # --- maximum of pseudo-likelihood function
  mlpl.phi <- phi.mnlogit$logLik

  # --- Information criteria for pseudo-likelihood
  nparm <- nless*nitems +  Maxnphi
  AIC <- -1*mlpl.phi - nparm
  BIC <- -2*mlpl.phi - nparm*log(length(unique(Master$PersonID)))

  results <- list(estimates=estimates,
                  fstack = fstack,
                  phi.mnlogit = phi.mnlogit,
                  mlpl.phi= mlpl.phi,
                  estimates = estimates,
                  Phi.mat = Phi.mat,
                  AIC = AIC[1],
                  BIC = BIC[1]
  )

  return(results)
}
