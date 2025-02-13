% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitStack.R
\name{FitStack}
\alias{FitStack}
\title{Up-dates association parameters of the nominal model}
\usage{
FitStack(
  Master,
  item.log,
  phi.log,
  fstack,
  TraitByTrait,
  pq.mat,
  npersons,
  nitems,
  ncat,
  nless,
  ntraits,
  Maxnphi,
  PhiNames,
  LambdaNames
)
}
\arguments{
\item{Master}{Master data set from which stacked data is created}

\item{item.log}{Last row contains current scale values (item.history)}

\item{phi.log}{Last row contains current estimates of phi}

\item{fstack}{Formula for stacked regression}

\item{TraitByTrait}{inTraitAdj matrix}

\item{pq.mat}{Summing array to get rest scores and totals}

\item{npersons}{Number of persons}

\item{nitems}{Number of items}

\item{ncat}{Number of categories per item}

\item{nless}{Number of categories less 1 (unique lambdas & unique nus)}

\item{ntraits}{Number of latent traits}

\item{Maxnphi}{Number of phis to be estimated}

\item{PhiNames}{Names of the Phi parameters}

\item{LambdaNames}{Names of lambdas that correspond to those in Master}
}
\value{
Phi.mat		 Matrix of up-dated estimates of assocation (phi) parameters

phi.log 	 History of iterations log likelihood and estimates of
                      lambda and phi parameters
}
\description{
Discrete choice model (conditional multinominal logistic regression model)
is fit to stacked data to up-date the matrix of association parameters of
the LMA that corresponds to the nominal item response model. This is a
function internal to 'fit.nominal' and is used for multi-dimensional models.
The function is similar to 'fit.StackGPCM'. This function is unlikely to
be run outside of 'fit.nominal' or 'ple.lma'.
}
