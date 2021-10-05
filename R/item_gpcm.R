#' Estimates item parameters of LMA with linear restrictions on category scores
#'
#' This function is internal to the function 'fit.gpcm' and performs the item
#' regressions. It is a core function of the pseudo-likelihood algorithm for
#' items of the GPCM. The function calls function 'itemGPCM.data' to create
#' the data for input into 'mlogit', which is use to fit a conditional
#' multinomial model for each item.  The up-dated scale values are put into
#' the Master data frame and the 'item.log' array.  It generally would not
#' run outside of 'fit.gpcm' or 'ple.lma'.
#'
#' @param	Master			 Master data frame
#' @param item.log		 History over iterations of items' log likelihood and
#'                            estimates of lambda, and item parameters
#' @param	Phi.mat			 Starting value of matrix of association parameters (optional)
#' @param TraitByTrait Trait adjacency matrix (same as inTraitAdj)
#' @param fitem		    	Formula for item regressions
#' @param PersonByItem  Same as inData
#' @param npersons	  	Number of persons
#' @param nitems			  Number of items
#' @param ncat			    Number of categories per item
#' @param nless			    Number of unique lambdas and unique nus per item
#' @param ntraits			  Number of latent traits
#' @param Maxnphi			  Number of phi parameters to bet estimated (NULL for 1 dimensional)
#' @param pq.mat			  Used to compute rest-scores and totals
#' @param	starting.sv 	Fixed category scores
#' @param LambdaName 		Lambda names for formula for items item regressions
#'
#' @return	Master     	Master data frame with up-dated category scores for items
#' @return	item.log   	Up-dated history array over iterations of the algorithm
#'                        of items' log likelihood and estimated lambda and
#'                        alpha parameters
#'
#' @export
item.gpcm <- function(Master, item.log, Phi.mat, fitem, TraitByTrait,
                      PersonByItem,  npersons, nitems, ncat, nless,
                      ntraits, Maxnphi, pq.mat,
                      starting.sv, LambdaName) {

 for (item in 1:nitems) {
    DataForItem  <-  ItemGPCM.data(Master, ItemID=item, Phi.mat, TraitByTrait,
                                  pq.mat, starting.sv, npersons, nitems, ncat,
                                  nless,  ntraits, LambdaName)

    data.fit  <- DataForItem$ItemFit
    xij       <- DataForItem$xij
    gpcmitem  <-  dfidx::dfidx(data.fit, choice="y", idx=c("CaseID","alt"))
    model.fit <- mlogit::mlogit(fitem, gpcmitem)

# --- Save output
    parms <- model.fit$coefficients[1:(nless+1)]    
    logLike <- as.numeric(model.fit$logLik)
    item.log[[item]] <- rbind(item.log[[item]],c(item,logLike,parms))

# --- Faster than my first method
   a <- parms[length(parms)]
   p1 <- 1
   p2 <- nitems*ncat
   for (person in 1:npersons) {
       this.column <-  5 + nitems*nless + item
       Master[p1:p2,this.column] <- a * xij[PersonByItem[person,item]]
       p1 <- p2 + 1
       p2 <- p1  + nitems*ncat - 1
   }

 }
	itemGPCM.results <- list(Master = Master, item.log = item.log)
	return(itemGPCM.results)
}
