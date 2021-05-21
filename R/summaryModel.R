#' Prints summary information about nominal or GPCM fit to data
#'
#' The output from fitting a model contains lots of information as a
#' list.  This function take model fit object and returns a summary
#' of model fit that is printed in console.  Information includes
#' model specification, information about data, convergence information,
#' and global fit statistics. This function is designed for Nominal,
#' generalized partial credit, and Rasch models.
#'
#' @param   opleLMA Object with results from ple.lma for "nominal", "gpcm", or "rasch"
#'
#'
#' @export

summaryModel <- function(opleLMA) {
      model.type  <- opleLMA[[1]]
      ntraits     <- opleLMA[[14]]
      nitems      <- opleLMA[[10]]
      ItemByTrait <- opleLMA[[3]]
	    TraitByTrait<- opleLMA[[2]]
	    ncat        <- opleLMA[[11]]
	    nless       <- opleLMA[[12]]
	    ItemNames   <- opleLMA[[5]]
	    npersons    <- opleLMA[[9]]
			tol         <- opleLMA[[16]]
			criterion   <- opleLMA[[17]]
			Maxnphi     <- opleLMA[[13]]
			mlpl.item   <- opleLMA[[24]]
			item.log    <- opleLMA[[18]]


cat("Summary of Output from pleLMA (beta version)   \n \n")

cat("Model and Data Information: \n")
	cat("  Model: ", model.type, "\n")
	cat("  Number of dimensions: ", ntraits, "\n")
	cat("  Number of persons:", npersons, "\n")
	cat("  Number of items: ", nitems, "\n")
	cat("  Number of categories per item:", ncat, "\n")
	if (model.type=="nominal" ) {
	   n.parms <- 2*nless*nitems + Maxnphi - ntraits
    }
  if (model.type=="gpcm") {
     n.parms <- nless*nitems + nitems + Maxnphi - ntraits
  }
  if (model.type=="rasch") {
     n.parms <- nless*nitems + Maxnphi
  }

  cat ("  Total number of unique parameters estimated:", n.parms, "\n")
	cat("     Number of unique marginal effects per item:", nless,"\n ")
	if (model.type=="nominal") {
	cat("     Number of unique nus per item:" , nless, "\n")
	}
  if (model.type=="rasch") {
  cat("     Number of unique assocation parameters", Maxnphi, "\n \n")

  } else {
	cat("     Number of unique assocation parameters", Maxnphi-ntraits, "\n \n")
  }
	cat("  Trait x Trait Adj Matrix: \n")
	print(TraitByTrait)
	cat("\n")

	cat("  Item x Trait Adj Matrix: \n")
	rownames(ItemByTrait) <- ItemNames
	print(ItemByTrait)
    cat("\n")


cat("Convergence: \n")
 if (model.type == "rasch") {
   n.iter <- opleLMA$phi.mnlogit$est.stat$niters
   criterion <- opleLMA$phi.mnlogit$est.stat$funcDiff
   tol <-  opleLMA$phi.mnlogit$est.stat$stopCond
   mlpl.item <- opleLMA$phi.mnlogit$logLik

 } else {
  	hist <-item.log[[1]]
    n.iter <- nrow(hist)-2
 }
  	cat("  Number of iterations: ", n.iter  , "\n")
    cat("  Tolerance Succesive loglik difference < ftol (", tol,") \n")
  	cat("  Criterion", criterion,"\n")
  	cat("\n")

cat("Model fit Statistics: \n")
    cat("  Maximum Log Pseudo-likelihood:", mlpl.item,"\n")
  	cat("  AIC:  ",(-2*mlpl.item + n.parms),"\n")
  	cat("  BIC:  ",(-2*mlpl.item + n.parms*log(npersons)),"\n")

}

