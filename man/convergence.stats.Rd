% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convergence.stats.R
\name{convergence.stats}
\alias{convergence.stats}
\title{Computes statistics to assess convergence of the nominal model}
\usage{
convergence.stats(item.log, nitems, nless, LambdaName, NuName)
}
\arguments{
\item{item.log}{Iteration history of items' log likelihoods and parameter estimates}

\item{nitems}{Number of items}

\item{nless}{Number of unique marginal terms (i.e., lambdas) and
unique category scale values(i.e., nus)}

\item{LambdaName}{Names of lambdas in item regressions}

\item{NuName}{Names of nu in item regressions}
}
\value{
diff.last   	 Differences between item loglikes & parameters on last two iterations

criterion.loglike  Maximum over items of the absolute value of LogLike differences

criterion.items  Sum of item differences of item parameters
}
\description{
For the nominal model, convergence statistics are computed for each item, as
well as the algorithm as a whole. The main argument is the history or log
from fitting item regressions. The convergence statistics are the differences
between current values of the log likelihoods and item parameter estimates
and those from the previous iteration. The maximum over item of the
differences of the log likelihood values is used to determine convergence
of the pseudo-likelihood algorithm.  This function is used internally, but
it can also be used after fitting a model to examine how many iterations
are required for parameter estimates to get close to the final values
and whether any item parameters are still changing.
}
\examples{
\donttest{
 # 9 items from dass data for 250 cases
 data(dass)
 inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]

#---   input for uni-dimensional
 inTraitAdj  <- matrix(1, nrow=1, ncol=1)
 inItemTraitAdj <- matrix(1, nrow=9, ncol=1)

 #--- Uni-dimensional Nominal Model
 n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj,inTraitAdj, tol=1e-02)

# Since this function in internal to fit.nominal, need to also run
 s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj)

 convergence.stats(n1$item.log, n1$nitems, n1$nless, s$LambdaName, s$NuName)

#--- Multidimensional models
#--- re-define inTraitAdj and inItemTraitAdj for 3 dimensional models
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)

  dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
  anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
  stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
  das <- list(dpress, anxiety, stress)
  inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])

#--- 3 dimensional nominal
  n3 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
  s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj)
  convergence.stats(n3$item.log, n3$nitems, n3$nless, s$LambdaName, s$NuName)
}

}
