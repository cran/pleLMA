#' loops through items and up-dates estimates of scale values for each item in Nominal Model
#'
#' This is a core function of the pseudo-likelihood algorithm for items of
#' the nominal model. The function calls function 'ItemData' to create
#' the data frame for input into 'mlogit', which is use to fit a conditional
#' multinomial model (i.e., a discrete choice model) for each item.  The
#' up-dated scale are put into the Master data frame and added to the
#' item.log array.  Generally the function would not run outside of
#' 'fit.nominal' or 'ple.lma'.
#'
#' @param 	Master			Current master frame
#' @param 	item.log	    Iteration history for the items parameters
#' @param   Phi.mat			Current estimate of Phi.mat
#' @param   PersonByItem    Person by item adjacency matrix (same as inData)
#' @param   npersons        Number of persons
#' @param   nitems          Number of items
#' @param   ntraits         Number of traits
#' @param   ncat            Number of categories
#' @param   nless           Number of unique lambda and unique nus (ncat-1)
#' @param   TraitByTrait    TraitsByTrait adjacency matrix (sam as TraitAdj)
#' @param   pq.mat          One dimensional array for computing rest-scores
#' @param   LambdaName      Marginal effect names used in formula and item
#'                             data frame for item regressions
#' @param   NuName          Scale values names in used in formula and item
#'                             data frame for item regressions
#' @param   fitem           Formula for item regression
#'
#'
#' @return	Master	  Master data frame up-dated scale values for all items
#' @return	item.log  Iteration history of item parameters where the last
#'                        row showing results from the current iteration
#'
#' @export
##################################################################################
ItemLoop <- function(Master, item.log, Phi.mat=Phi.mat, PersonByItem, npersons,
                     nitems, ntraits, ncat, nless, TraitByTrait, pq.mat,
                     LambdaName, NuName, fitem) {
for (item in 1:nitems) {
  	DataForItem  <- ItemData(Master, ItemID=item, Phi.mat=Phi.mat, npersons,
  		                         nitems, ntraits, ncat, nless, TraitByTrait,
  		                         pq.mat, LambdaName, NuName)

	nomitem  <-  dfidx::dfidx(DataForItem, choice="y", idx=c("CaseID","alt"))

   	model.fit    <- mlogit::mlogit(fitem, nomitem)

# --- save outpupt to log
    parms <- model.fit$coefficients[1:(2*nless)]

# --- maximum of log pseudo-likelihood function
    mlpl.item <- as.numeric(model.fit$logLik)

  	item.log[[item]] <- rbind(item.log[[item]],
	                          c(item, mlpl.item, parms))

# --- Faster than my first method
   	tmp <- parms[(ncat):(2*nless)]
   	NuEst <- c(-sum(tmp),tmp)       # compute nu_1
   	p1 <- 1
   	p2 <- nitems*ncat
  	for (person in 1:npersons) {
    	this.column <-  5 + nitems*nless + item
      	Master[p1:p2,this.column] <- NuEst[PersonByItem[person,item]]
   		p1 <- p2 + 1
   		p2 <- p1  + nitems*ncat - 1
   		}
}
ItemLoop.Results <- list(Master=Master,item.log=item.log)
return(ItemLoop.Results)
}
