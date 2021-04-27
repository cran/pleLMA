#' Imposes scaling constraint to identify parameters of LMA (GPCM)
#'
#' Scaling is internal to the function 'fit.gpcm', which fits the GPCM
#' version of the LMA. It imposes the required scaling identification
#' constraint by transforming the conditional covariance matrix 'Phi.mat' to a
#' conditional correlation matrix. The inverse transformation is applied to the
#' current estimates of the slope or 'a' parameters. Category scale values are
#' recomputed using the re-scale slopes (i.e., nu= a*x) and these are put back
#' into the Master data set so that data are ready for the next iteration of the
#' algorithm.
#'
#' @param   Master	 	    Master/main data set
#' @param   item.log      Iteration history array, last row are current parameters
#' @param   Phi.mat		    Current phi matrix
#' @param   PersonByItem 	inData (response patterns)
#' @param   npersons      Number of persons
#' @param  	nitems		    Number of items
#' @param   ncat 		    	Number of categories
#' @param   nless			    Number of unique lambdas (ncat-1)
#' @param   ntraits	    	Number of latent traits
#' @param   starting.sv		Matrix of fixed category scores (nitems x ncat)
#' @param   item.by.trait Object that indicates which trait item loads on
#'
#' @return  Master   Master data set with re-scaled scale values
#' @return  Phi.mat	 Re-scaled matrix of association parameters
#'
#' @export

ScaleGPCM <- function(Master, item.log, Phi.mat, PersonByItem, npersons, nitems,
                      ncat, nless, ntraits, starting.sv, item.by.trait) {

# --- Re-scale phi so that it's a correlation matrix (i.e., set p[p,p]=1)

    idPhi <- solve( diag(diag(Phi.mat))  )
    c <- sqrt(idPhi)
    Phi.mat <- c %*% Phi.mat %*% c          # This is new Phi

# --- get a from last line of item.log
    All.a <- matrix(NA,nrow=nitems,ncol=ncat)
    for (item in 1:nitems) {
      tmp <- item.log[[item]]
      All.a[item,] <- matrix(tmp[nrow(tmp), ncol(tmp)], nrow=1, ncol=1)
    }

# --- Rescale the alphas  (could be more efficient but gets job done)
    sNu.p <- matrix(0,nrow=ncat,ncol=nitems)
    for (p in 1:ntraits) {
      for (i in 1:nitems) {
       if (item.by.trait[,i] == p) {
         sa   <- All.a[i,]/c[p,p]
         sNu.p[ , i] <- as.vector(sa) * as.matrix(starting.sv[i, ])
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

   results<- list(Master, Phi.mat=Phi.mat)
   return(results)
}
