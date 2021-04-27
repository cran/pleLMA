#' Prepares data for up-dating scale value parameters of nominal model
#'
#' This function creates a data frame, 'item data', to be used in the item
#' regressions for nominal models. It computes weighted rest scores and totals,
#' including correlated traits. This function is internal to 'ItemLoop' and
#' it is unlikely to be run outside of 'fit.nominal' or 'ple.lma'.
#'
#' @param  	Master  	   Master data frame
#' @param  	ItemID		   The item for which scale values are being up-dated
#' @param 	Phi.mat		   Current estimate conditional covariance matrix
#'                           (i.e., association paramters)
#' @param   npersons 	   Number of persons
#' @param   nitems		   Number of items
#' @param   ntraits		   Number of traits
#' @param   ncat		     Number of categories
#' @param   nless		     Number of unique lambdas and unique nus
#' @param   TraitByTrait Same as inTraitAdj
#' @param   pq.mat	     One dimemsinal array used to get rest and totals scores
#' @param   LambdaName	 Name of lambdas for in item regression
#' @param   NuName	     Name of nus in item regression
#'
#' @return	ItemFit 	Data frame used to up-date scale values
#'
#' @export
ItemData <- function(Master, ItemID, Phi.mat=Phi.mat, npersons, nitems,
                     ntraits, ncat, nless, TraitByTrait, pq.mat, LambdaName,
                     NuName) {
	item.data<- Master[which(Master$Item==ItemID),]

  # --- Do lambda.select only once, make global and put them in matrix
    	jcol <- 1
     	lambda.select <- matrix(NA,nrow=1,ncol=nless)
    	for (item in 1:nitems) {
	    	if (item == ItemID) {
     		  for (cat in 1:nless) {
        	   lambda.select[1,(jcol)] <- paste("lam", item, 'j', (cat+1), sep="")
         		jcol <- jcol + 1
       		}
	 	  }
   	}
	effects <- as.matrix(item.data[ , which(names(item.data) %in% lambda.select)])

  # --- Compute Rest & test totals for item data
 	nu.ItemID <- as.matrix(item.data[, (5+nitems*nless+1):(ncol(item.data)-2)])
	wrest <- matrix(0, nrow=ncat*npersons, ncol=1)
	trait.combo <- 1
    	for (p in 1:ntraits) {
	      for (q in p:ntraits) {
	        if (TraitByTrait[p,q] == 1) {
	          pq <- as.matrix(pq.mat[ ,ItemID, trait.combo] )
       	  	  wrest <- wrest + (nu.ItemID  %*% pq)  * Phi.mat[p, q]
	      	  trait.combo <- trait.combo + 1
	        }
	      }
    	}
	Nu <- matrix(NA,nrow=(npersons*ncat),ncol=nless)
	for (cat in 1:nless) {
   	Nu[,cat] <-  wrest * effects[,cat]
	}

# --- Build matrix to input to mnlogit package
	ItemFit <- cbind(item.data[, c(2,4,5)], effects, Nu, item.data$alt,
	                 item.data$choice)

# ---data used to fit model to data
	ItemFit <- as.data.frame(ItemFit)
    names(ItemFit) <- c("CaseID", "Category", "y", LambdaName, NuName,"alt" ,
                        "choice")

    return(ItemFit)
}
