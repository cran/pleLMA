#' Produces a summary of results
#'
#' This utility function creates a summary list with five elements.
#' The first is a 'report' that contains a summary of characteristics
#' of the data, the model specification, convergence information,
#' and fit statistics. The second and third elements complete the
#' model specification and are the trait adjacency matrix
#' ('TraitByTrait') and  the item x trait adjacency matrix
#' ('ItemByTrait'), respectively. The fouth element, 'estimates',
#' contains a matrix of item parameters, and the fifth element,
#' 'phi.mat' contains association parameter estimates.
#'
#' @param   model.fit  A list object from fitting a model to data
#'
#' @return  results   A list with summary information
#'
#' @examples
#' #--- 3 items from depression, anxiety and stress scales of
#' #    the daass for 250 cases
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#' #---  log-linear model of independence
#'  ind <- ple.lma(inData, model.type="independence")
#'  noquote(lma.summary(ind))
#'
#' #---   input for uni-dimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#' #--- rasch family
#'  r1 <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
#'  lma.summary(r1)
#'  #--- Or if specific output is desired
#'  noquote(lma.summary(r1)$report)
#'  lma.summary(r1)$TraitByTrait
#'  lma.summary(r1)$ItemByTrait
#'  lma.summary(r1)$estimates
#'  lma.summary(r1)$phi
#' 
#' \donttest{
#' #--- generalized parial credit model
#' g1 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1e-03)
#' lma.summary(g1)$report
#' lma.summary(g1)$estimates
#' lma.summary(g1)$phi
#'
#' #--- nominal response model
#' n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
#' noquote(lma.summary(n1))
#' }
#'
#' @export
lma.summary <- function(model.fit) {

 #--- read in needed output from model
  model.type  <- model.fit$model.type
  TraitByTrait<- model.fit$TraitByTrait
  ItemByTrait <- model.fit$ItemByTrait
  ItemNames   <- model.fit$ItemNames
  formula.item<- model.fit$formula.item
  formula.phi <- model.fit$formula.phi
  npersons    <- model.fit$npersons
  nitems      <- model.fit$nitems
  ncat        <- model.fit$ncat
  nless       <- model.fit$nless
  Maxnphi     <- model.fit$Maxnphi
  ntraits     <- model.fit$ntraits
  starting.sv <- model.fit$starting.sv
  tol         <- model.fit$tol
  criterion   <- model.fit$criterion
  item.log    <- model.fit$item.log
  phi.log     <- model.fit$phi.log
  estimates   <- model.fit$estimates
  Phi.mat     <- model.fit$Phi.mat
  item.mnlogit<- model.fit$item.mnlogit
  phi.mnlogit <- model.fit$phi.mnlogit
  mlpl.item   <- model.fit$mlpl.item
  mlpl.phi    <- model.fit$mlpl.phi
  AIC         <- model.fit$AIC
  BIC         <- model.fit$BIC

# Model specific information

  if (model.type=="independence") {
         mlpl <- mlpl.phi
   }
  if (model.type=="rasch") {
         mlpl <- mlpl.phi
  }
  if (model.type=="gpcm") {
         mlpl <- mlpl.item
  }
  if (model.type=="nominal") {
         mlpl <- mlpl.item
  }

  if (model.type == "rasch" || model.type=="independence") {
     n.iter <- 1
     criterion <- phi.mnlogit$est.stat$funcDiff
     tol <-  phi.mnlogit$est.stat$stopCond

} else {
     hist <-item.log[[1]]
     n.iter <- nrow(hist)-2
}

   n.lambda.unique <- nless*nitems

  if (model.type=="nominal" ) {
      n.parms <- 2*nless*nitems + Maxnphi - ntraits
      n.nu <- nless*nitems
      n.phi.unique <- Maxnphi-ntraits
  } else if (model.type=="gpcm") {
      n.parms <- nless*nitems + nitems + Maxnphi - ntraits
      n.nu <- nitems
      n.phi.unique <- Maxnphi-ntraits
  } else if (model.type=="rasch") {
      n.parms <- nless*nitems + Maxnphi
      n.nu <- 0
      n.phi.unique <- Maxnphi
  } else if (model.type=="independence") {
      n.parms <- nless*nitems
      n.nu <- 0
      ntraits <- 0
      Maxnphi <- 0
      n.phi.unique <- 0
      Phi.mat <- NULL
      TraitByTrait <- NULL
      ItemByTrait <- NULL

  }



date_time <- Sys.time()

  report <- rbind(
    "                                                         ",
    "=========================================================",
    paste("Pseudo-likelihood Estimation of", model.type, "model"),
    "=========================================================",
    paste("Report Date: ", date_time ),
    "                                                         ",
    paste("Data Information:"),
    paste("   Number of cases/individuals ", npersons),
    paste("   Number of items ", nitems),
    paste("   Number of categories per item", ncat),
    paste("   Number of dimensions: ", ntraits),
    "                                                        ",
    paste("Model Specification:"),
    paste("  Number of unique parameters", n.parms),
    paste("  Number of unique marginal effects: ",  n.lambda.unique),
    paste("  Number of unique category parameters (nu's or a's): ", n.nu),
    paste("  Number of unique association parameters (phis):",  n.phi.unique),
    "                                                         ",
    paste("Convergence Information:"),
    paste("  Number of iterations: ", n.iter),
    paste("  Tolerence set tol", tol),
    paste("  Criterion ", criterion),
    "                                                         ",
    paste("Model Fit Statistics:"),
    paste("  Maximum log pseudo-likelihood function:", mlpl),
    paste("  AIC:  ", AIC),
    paste("  BIC:  ", BIC),
    "                                                         "
    )


   report.all <- list(report = report,
                      TraitByTrait = TraitByTrait,
                      ItemByTrait = ItemByTrait,
                      estimates = estimates,
                      phi = Phi.mat)
   results <- (report.all)
  return(results)
}





