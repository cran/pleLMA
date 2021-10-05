#' Creates data frame up-dating phi parameters of the gpcm.
#'
#' This function creates a data frame, 'gpcm.item.data', to be used in the item
#' regressions of LMA models where category scales values are fixed. Sets up data
#' for up-dating alpha parameters of the LMA that corresponds to the GPCM. This
#' function is internal to 'item.gpcm' and it is unlikely to be run outside of
#' 'fit.gpcm' or 'ple.lma'.
#'
#' @param Master	Data frame of all data in long format
#' @param ItemID  Specifies the item for which a data frame is created to be
#'                 input into an item regression
#' @param Phi.mat	 Starting value of matrix of association parameters
#' @param TraitByTrait Trait by trait adjacency matrix (same as inTraitAdj)
#' @param pq.mat	Array used to compute rest scores
#' @param starting.sv   Matrix of item category scores that are fixed
#' @param npersons	Number of persons
#' @param nitems	Number of items
#' @param ncat		Number of categories
#' @param nless		Number of unique lambdas (ncat-1)
#' @param ntraits Number of latent traits
#' @param LambdaName  Names for lambda for item regression
#'
#' @return gpcm.item.data  Data frame for item to be used up-dated in an item regression for specified item
#'
#' @export

ItemGPCM.data <- function(Master, ItemID, Phi.mat, TraitByTrait, pq.mat, starting.sv,
                          npersons, nitems, ncat, nless, ntraits, LambdaName) {

    item.data<- Master[which(Master$Item==ItemID),]

# --- Do lambda.select only once, make global and put them in matrix
    jcol <- 1
    lambda.select <- matrix("name",nrow=1,ncol=nless)
    for (item in 1:nitems) {
	  if (item == ItemID) {
       for (cat in 1:nless) {
         lambda.select[1,(jcol)] <- paste("lam",item,'j',(cat+1), sep="")
         jcol <- jcol + 1
       }
	  }
   }
	effects <- as.matrix(item.data[ , which(names(item.data) %in% lambda.select)])

# --- Rest & test totals for item data to estimate new alpha
	nu.ItemID <- as.matrix(item.data[, (5+nitems*nless+1):(ncol(item.data)-2)])
	alpha <- matrix(0, nrow=ncat*npersons, ncol=1)
	trait.combo <- 1
    for (p in 1:ntraits) {
	   for (q in p:ntraits) {
	     if (TraitByTrait[p,q] == 1) {
	       tmp <- pq.mat[ , ItemID, trait.combo]
	       pq <- as.matrix(tmp)
         alpha <- alpha + (nu.ItemID  %*% pq)  * Phi.mat[p, q]
		     trait.combo <- trait.combo + 1
	     }
	   }
    }

# --- the weights (could be done once but will repeat for now)
   	xij <- matrix(rep(starting.sv[item,], npersons), nrow=ncat*npersons, ncol=1)
    alpha <- alpha * xij

# --- Build matrix to input to mnlogit package
    ItemFit <- cbind(item.data[, c(2,4,5)], effects, alpha, item.data$alt, item.data$choice)

# ---data used to fit model to data
    ItemFit <- as.data.frame(ItemFit)
    names(ItemFit) <- c("CaseID", "Category", "y", LambdaName, "alpha","alt","choice")

    gpcm.item.data <- list (ItemFit= ItemFit,
                            xij = xij)
    return(gpcm.item.data)
}
