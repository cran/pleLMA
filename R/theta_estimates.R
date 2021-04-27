#' Computes estimates of theta (values on latent trait(s)) for all LMA models
#'
#' The final estimates of the item scale values and the conditional
#' covariance matrix (i.e, Phi.mat) are used to compute values on latent 
#' traits for each individual or case.  The estimated thetas are the 
#' (conditinal) mean values of response patterns. The correlations 
#' between the estimated thetas equal the marginal correlations. 
#'
#' @param  model.fit     Object containing output from running ple.lma
#' @param  inData        Matrix of response patterns
#'
#' @return	theta.est 	 A person by trait matrix of values on the latent traits
#'
#' @examples
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#'  r1 <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
#'  theta.r1 <- theta.estimates(inData, r1)
#'
#' \donttest{
#' g1 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)
#' theta.g1 <- theta.estimates(inData, g1)
#'
#' n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj,inTraitAdj)
#' theta.n1 <- theta.estimates(inData, n1)
#' }
#'
#' @export
theta.estimates <- function(inData,model.fit) {

  # --- pull information needed
  PersonByItem <- inData
  estimates <- model.fit$estimates
  ItemByTrait <- model.fit$ItemByTrait
  model.type <- model.fit$model.type
  Phi.mat <- model.fit$Phi.mat
  npersons <- model.fit$npersons
  nitems <- model.fit$nitems
  ncat <- model.fit$ncat
  ntraits <- model.fit$ntraits
  phi.mnlogit <- model.fit$phi.mnlogit

  TraitByItem <- t(ItemByTrait)
  theta.est <- matrix(0,nrow=npersons,ncol=ntraits)
  trait.by.nu <- matrix(NA, nrow=ntraits, ncol=nitems)

  if (model.type=="nominal") {
    for (person in 1:npersons) {
      nu.est <- estimates[,(2+ncat):ncol(estimates)]
      for (trait in 1:ntraits) {
        for (item in 1:nitems) {
          trait.by.nu[trait,item] <-
            TraitByItem[trait, item] * nu.est[item, PersonByItem[person,item]]
        }
      }
      theta.est[person, ] <- rowSums ( Phi.mat %*% trait.by.nu )
    }
  }

  if (model.type=="gpcm") {
    a.est <- estimates[,(2+ncat)]
    nu.est <- a.est * estimates[, (3+ncat):ncol(estimates)]
    for (person in 1:npersons) {
      for (trait in 1:ntraits) {
        for (item in 1:nitems) {
          trait.by.nu[trait,item] <-
            TraitByItem[trait, item] * nu.est[item, PersonByItem[person,item]]
        }
      }
      theta.est[person, ] <- rowSums ( Phi.mat %*% trait.by.nu )
    }
  }

  if (model.type=="rasch") {
    nu.est <- estimates[,(ncat+1):ncol(estimates)]
    phi <- Phi.mat
    pointer <- 1
    phi.mtx <- matrix(0, nrow=ntraits, ncol=ntraits)
    for (p in 1:ntraits){
      for (q in p:ntraits) {
        if (p==q) {
          phi.mtx[p,p] <- phi[pointer]
          pointer <- pointer + 1

        } else {
          phi.mtx[p,q] <- phi[pointer]
          phi.mtx[q,p] <- phi[pointer]
          pointer <- pointer + 1
        }
      }
    }
    for (person in 1:npersons) {
      for (trait in 1:ntraits) {
         for (item in 1:nitems) {
             trait.by.nu[trait,item] <-
              TraitByItem[trait, item] * nu.est[item, PersonByItem[person,item]]
      }
    }
    theta.est[person, ] <- rowSums ( phi.mtx %*% trait.by.nu )
  }
  }
  return(theta.est = theta.est)

}

