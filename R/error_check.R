#' Checks for basic errors in input to the 'ple.lma' function
#'
#' This functions looks at the input to the main function (ple.lma) and checks
#' for 11 different possible errors. If an error is detected, the function
#' issues a warning and stops any further execution.  This funcion is internal
#' to 'ple.lma' but can be used outside of the wrapper function.
#'
#' @param inData         Data frame with columns corresponding to categorical variables
#'                       and rows to the number of cases
#' @param model.type     Type of model that will be fit to data
#' @param inTraitAdj     Trait x Trait adjacency matrix (not required for independence)
#' @param inItemTraitAdj Item x Trait adjacency matrix (not required for independence)
#'
#' @return Message whether error was detected in input, and if so the nature of the error
#'
#' @examples
#'  #--- some data
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#'  #--- no errors
#'  error.check(inData, model.type="independence")
#'
#'  #--- for unidimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#'  #--- no errors
#'  error.check(inData, model.type="rasch", inTraitAdj, inItemTraitAdj)
#'  error.check(inData, model.type="gpcm", inTraitAdj, inItemTraitAdj)
#'  error.check(inData, model.type="nominal", inTraitAdj, inItemTraitAdj)
#'
#'
#' @export
######################################################################################
error.check <- function(inData, model.type, inTraitAdj=NULL, inItemTraitAdj=NULL) {

if (dim(table(is.na(inData))) > 1) {
   stop("Your data may contain missing values.  This code is not set up to handle missing values.")
	}

if (min(inData) != 1) {
    stop("The responses to items should be consecutive integers that run from 1 to number of categories")
}


"%notin%" <- Negate("%in%")
model.types <- c("independence", "rasch", "gpcm", "nominal")
if ( model.type %notin% model.types) {
    stop("Model type not recognized.  Possible types are 'independence', 'rasch', 'gpcm' and 'nominal' ")
}

#--- Errors specific to rasch, gpcm and nominal models

if (model.type != "independence") {

  if (length(inTraitAdj) == 0) {
    stop("Missing inTraitAdj")
  }

  if (length(inItemTraitAdj) == 0) {
    stop("Missing inItemTraitAdj")
  }

  for (i in 1:nrow(inItemTraitAdj)) {
    numbers <- sort(unique(inData[,i]))
    diffs <- diff(numbers)
    if (max(diffs) >1 ) {
     stop("The responses in your data set are not consecutive integers")
  }
}


  if (nrow(inTraitAdj) != ncol(inTraitAdj)) {
    stop("The Trait x Trait adjacency matrix must be square")
	}

  if (nrow(inTraitAdj) != ncol(inItemTraitAdj)) {
   stop("The Item x Trait matrix must have same number of columns as the Trait x Trait Matrix")
   }

  if (max(rowSums(inItemTraitAdj)) > 1) {
   stop("This code is only set up for simple structure; that is, each item can only load on one item.  You should check you Item x Trait matrix")
   }

  if (sum(inItemTraitAdj) != nrow(inItemTraitAdj)) {
  stop("Check your Item x Trait matrix.  There should be one 1 in each row and other elements should be 0")
   }


Two <- c("gpcm","nominal")
if (nrow(inItemTraitAdj)==2) {
  if (model.type %in% Two) {
    stop("You only have have 2 variables, consider using package logmult or gnm. Currently only model types independence and rasch can handle two variables.")
  }
 }
}
print("No errors detected in the input")
}
