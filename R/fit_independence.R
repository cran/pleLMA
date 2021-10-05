#' Fits the log-linear model of independence  
#'
#' This function fits by the log-linear model of independence (i.e., only
#' includes marginal effect terms) using pseudo-likelihood estimation.
#' This provides a baseline model with which to compare other models.
#' The independence maximumn of the loglikehood can be used is a measure
#' of no association. The input to the function is only the Master data
#' set and the names of marginal effect terms and items, all of which
#' are created by the 'set.up' function.  This function is called from
#' 'ple.lma' or can be run output of wrapper.
#'
#' @param  Master          Master data set from set.up
#' @param  LambdaNames     Needed to define formula
#' @param  LambdaName      Used for column names of matrix estimates
#' @param  ItemNames       Used for row names of number of item by parameter matrix of
#'                          estimated Lambda parameters
#'
#' @return phi.mlogit     Parameters estimates and mlpl = logLike output from mnlogit
#' @return fstack          Formual used in stacked regression
#' @return estimates       Item by parameter estimates matrix
#' @return mlpl.phi        Maximum of log pseudo-likelihood from stacked regression
#' @return AIC             Akaike information criterion for pseudo-likelihood (smaller is better)
#' @return BIC             Bayesian information criterion for pseudo-likelihood (smaller is better)
#'
#' @examples
#' #--- data and set-up
#' data(dass)
#' inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#' s <- set.up(inData, model.type='independence')
#'
#' #--- fit independence model
#' ind <- fit.independence(s$Master, s$LambdaNames, s$LambdaName, s$ItemNames)
#'
#' @export
fit.independence <- function(Master, LambdaNames, LambdaName, ItemNames) {

# --- make mlogit data
  master.mlogit <- dfidx::dfidx(Master, choice="y", idx=c("CaseID","alt"))

  fstack <- stats::as.formula(paste("choice ~ ", paste(LambdaNames,
            collapse = "+"), "| 0 | 0", sep=" ") )

  fit.stack <- mlogit::mlogit(fstack, master.mlogit)

 # --- save some things
    estimates <- matrix(fit.stack$coefficients,
	                 nrow=length(unique(Master$Item)),
					 ncol=(length(unique(Master$Category))-1),
					 byrow=TRUE)

    e1 <- -(rowSums(estimates))
    estimates <- cbind(e1,estimates)
    rownames(estimates) <- ItemNames
    colnames(estimates) <- c("lam1", LambdaName)

  # --- maximum of log pseudo-likelihood function
  mlpl.phi <- as.numeric(fit.stack$logLik)

  # --- AIC
  AIC <- -1*mlpl.phi - length(LambdaNames)

  # --- BIC
  BIC <- -2*mlpl.phi - length(LambdaNames)*log(length(unique(Master$PersonID)))

  summary.stack <- summary(fit.stack)

  results <- list(phi.mlogit = summary.stack,
                  fstack = fstack,
                  estimates = estimates,
                  mlpl.phi = mlpl.phi[1],
                  AIC = AIC[1],
                  BIC = BIC[1]
                  )

  return(results)
}
