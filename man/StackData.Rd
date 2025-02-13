% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StackData.R
\name{StackData}
\alias{StackData}
\title{Prepares data for up-dating association parameters of a multidimensional nominal LMA}
\usage{
StackData(
  Master,
  item.log,
  phi.log,
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
\item{Master}{Master data frame from which stacked data is created}

\item{item.log}{Last row contains current scale values (item.history)}

\item{phi.log}{Last row contains current estimates of phi}

\item{pq.mat}{Summing array to get rest scores and totals}

\item{npersons}{Number of persons}

\item{nitems}{Number of items}

\item{ncat}{Number of categories per item}

\item{nless}{Number of categories less 1 (i.e., unique lambdas and unique nus)}

\item{ntraits}{Number of latent traits}

\item{Maxnphi}{Number of phis to be estimated}

\item{PhiNames}{Names of the Phi parameters}

\item{LambdaNames}{Names of lambdas that correspond to those in Master}
}
\value{
Phi.mat		Up-dated matrix of phi parameters

phi.log 	History of iterations for LogLike, Lambda and phi parameters
}
\description{
Prepares data frame for input to 'mnlogit' for the stacked regression 
to obtain association parameters of the multidimensional LMA models
corresponding to the Nominal item response model. This function is
called from 'fit.nominal'. It generally would not run outside of 
either 'fit.nominal' or 'ple.lma'.
}
