#' Imposes scaling constraint to identify parameters LMA (nominal) model
#'
#' Scaling is internal to the function 'fit.nominal', which corresponds
#' the the nominal item response theory model.  It imposes the required
#' scaling identification constraint by transforming the conditional
#' covariance matrix 'Phi.mat' to a conditional correlation matrix (i.e.,
#' set phi_mm=1 for all m). The inverse transformation is applied to
#' the current category scale value estimates and these are put back
#' into the Master data frame so that data are ready for the next
#' iteration of the algorithm.
#'
#' @param   Master	      Current Master data frame.
#' @param   item.log      Iteration history of LogLike, lambda, and item parameters
#' @param   Phi.mat       Current phi matrix
#' @param   PersonByItem  inData
#' @param   npersons      Number of persons
#' @param   nitems        Number of items
#' @param   ncat          Number of categories
#' @param   nless         Number of unique nus (ncat-1)
#' @param   ntraits       Number of (latent) dimensions
#' @param   item.by.trait Indicates the trait an item load on.
#'
#' @return  Master       Master frame with re-scaled scale values
#' @return  Phi.mat		   Re-scaled matrix of association parameters
#'
#' @export
Scale <- function(Master, item.log, Phi.mat, PersonByItem, npersons, nitems, ncat,
                  nless, ntraits, item.by.trait) {

# --- Re-scale phi so that it's a correlation matrix (i.e., set p[p,p]=1)
    idPhi <- solve( diag(diag(Phi.mat))  )
    c <- sqrt(idPhi)
    Phi.mat <- c %*% Phi.mat %*% c          # This is new Phi

# --- get up-dated nus from last lines of item.history
    j1 <- 3 + nless
    j2 <- j1 + nless - 1
    AllNu <- matrix(0, nrow=nitems, ncol=nless)
    for (i in 1:nitems) {
      history<- item.log[[i]]
      AllNu[i,1:nless] <- history[nrow(history), j1:j2]
    }
    Nu1 <-  -rowSums(AllNu)
    Nu1 <- matrix(Nu1, nrow=nitems, ncol=1)
    AllNu <- cbind(Nu1,AllNu)

# --- Rescale the nus  (could be more efficient but gets job done
    sNu.p <- matrix(0,nrow=ncat,ncol=nitems)
    for (p in 1:ntraits) {
      for (i in 1:nitems) {
       if (item.by.trait[1,i] == p) {
         sNu.p[,i ] <- AllNu[i,]/c[p,p]
       }
     }
    }
# --- replace nus in Master with re-scaled ones
    d <- 5+nitems*nless
    p1 <- 1
    p2 <- nitems*ncat
    for (person in 1:npersons) {
       for (item in 1:nitems) {
         this.column <- d + item
         Master[p1:p2,this.column] <- sNu.p[PersonByItem[person,item],item]
       }
	 p1 <- p2 + 1
	 p2 <- p1  + nitems*ncat - 1
    }

   results<- list(Master, Phi.mat)
   return(results)
}
