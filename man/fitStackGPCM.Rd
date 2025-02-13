% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitStackGPCM.R
\name{fitStackGPCM}
\alias{fitStackGPCM}
\title{Up-dates association parameters of the GPCM by fitting model to stacked data}
\usage{
fitStackGPCM(
  Master,
  item.log,
  phi.log,
  fstack,
  TraitByTrait,
  starting.sv,
  npersons,
  nitems,
  ncat,
  nless,
  ntraits,
  Maxnphi,
  pq.mat,
  LambdaNames,
  PhiNames
)
}
\arguments{
\item{Master}{Master data set from which stacked data is created}

\item{item.log}{Needed to get most recent values of scale values (item.log)}

\item{phi.log}{History of estimates parameters from stacked regression}

\item{fstack}{Forumla for stacked regression}

\item{TraitByTrait}{inTraitAdj matrix}

\item{starting.sv}{Fixed category scores}

\item{npersons}{Number of persons}

\item{nitems}{Number of items}

\item{ncat}{Number of categories per item}

\item{nless}{Number of unique lambdas and unique nus per item}

\item{ntraits}{Number of latent traits}

\item{Maxnphi}{Number of phi parameters to bet estimated (NULL for 1 dimensional)}

\item{pq.mat}{Used to compute rest-scores and totals}

\item{LambdaNames}{Needed for formula and data for up-dating phi (stacked regresson)}

\item{PhiNames}{Null for 1D models}
}
\value{
Phi.mat	Up-dated matrix of phi parameters

item.log of iterations for LogLike, Lambda and phi parameters
}
\description{
Discrete choice model (conditional multinomial logistic regression) is fit to
stacked data to up-date matrix of association parameters of the LMA that
corresponds to the generalized partial credit model. This function is called
from 'fit.gpcm', which is called from 'ple.lma'.  It is unlikely that it
would be run outside of these wrappers.  It is only slightly different from
'fitStack' for nominal models.
}
