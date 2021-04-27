#'  Up-dates association parameters of the GPCM by fitting model to stacked data
#'
#' Discrete choice model (conditional multinomial logistic regression) is fit to
#' stacked data to up-date matrix of association parameters of the LMA that
#' corresponds to the generalized partial credit model. This function is called
#' from 'fit.gpcm', which is called from 'ple.lma'.  It is unlikely that it
#' would be run outside of these wrappers.  It is only slightly different from
#' 'fitStack' for nominal models.
#'
#' @param 	Master		Master data set from which stacked data is created
#' @param   item.log 	Needed to get most recent values of scale values (item.log)
#' @param   phi.log		History of estimates parameters from stacked regression
#' @param   fstack		Forumla for stacked regression
#' @param   TraitByTrait inTraitAdj matrix
#' @param   starting.sv  Fixed category scores
#' @param   npersons	Number of persons
#' @param   nitems		Number of items
#' @param   ncat	   	Number of categories per item
#' @param   nless		  Number of unique lambdas and unique nus per item
#' @param   ntraits		Number of latent traits
#' @param   Maxnphi		Number of phi parameters to bet estimated (NULL for 1 dimensional)
#' @param   pq.mat		Used to compute rest-scores and totals
#' @param   LambdaNames Needed for formula and data for up-dating phi (stacked regresson)
#' @param   PhiNames	Null for 1D models
#'
#' @return	Phi.mat	Up-dated matrix of phi parameters
#' @return	item.log of iterations for LogLike, Lambda and phi parameters
#'
#' @export

fitStackGPCM <- function(Master, item.log, phi.log, fstack, TraitByTrait,
                        starting.sv, npersons, nitems, ncat, nless, ntraits, Maxnphi,
                         pq.mat, LambdaNames, PhiNames) {

########--- First stack the regressions in one data set

  NewNu <- matrix(NA,nrow=nitems,ncol=ncat)
  for (item in 1:nitems) {
    tmp <- as.matrix(item.log[[item]])
    a <- matrix(tmp[nrow(tmp), ncol(tmp)], nrow=1, ncol=1)
    NewNu[item,] <- as.vector(a) * as.matrix(starting.sv[item,])
  }
  NewNu.block <- matrix(t(NewNu), nrow=nitems*ncat, ncol=1)    # New's for each cat & item

  # --- Repeat each block for each person,
  NuItemCat <- do.call(rbind, replicate(npersons, NewNu.block, simplify=FALSE))

  # --- Needed for matrices to conform
  NuItemCat.wide <- do.call(cbind, replicate(Maxnphi,NuItemCat,simplify=FALSE))

  # --- Need nperson x nu matrix of dim (npersons x nitems).
  #  Take 1st row from every nu block
  last.lambda <- 5+nitems*nless
  PersonByNu = as.matrix(Master[seq(1, nrow(Master), nitems*ncat), (last.lambda+1):(last.lambda+nitems)])

  # Compute values to estimate the phis
  phi.col <- matrix(0, nrow=nitems*npersons, ncol=Maxnphi)
  colnames(phi.col) <- PhiNames
  irow <- 1
  for (person in 1:npersons) {
    for (item in 1:nitems) {
      for (trait.combo in 1:Maxnphi) {
        phi.col[irow,trait.combo] <- as.vector(PersonByNu[person,] %*% pq.mat[ , item , trait.combo] )
      }
      irow <- irow + 1
    }
  }

  # --- Bump to Nstack
  Phi.col <- phi.col[rep(seq_len(nrow(phi.col)), each = ncat), ]

  # --- Weight rest.scores and totals by the nus for that item and category
  Phi.col <- as.data.frame(NuItemCat.wide * Phi.col)

  stack.data <- cbind(Master[,2:5], Master[,6:last.lambda], Phi.col, Master$alt, Master$choice)
  names(stack.data) <- c("CaseID", "Item", "Category", "y", LambdaNames, PhiNames, "alt", "choice")

  ####### ---Fit stacked data
	model.fit<- mnlogit::mnlogit(fstack, stack.data, choiceVar="Category")

  # --- Prepare output
	parms <- model.fit$coefficients
	logLike <- model.fit$logLik

	phi.log <- rbind(phi.log, c(logLike, parms))

	phis <- parms[(nless*nitems+1):(nless*nitems+Maxnphi)]
	Phi.mat <- matrix(0, nrow=ntraits, ncol=ntraits)
	iphi <- 1
	for (p in 1:ntraits) {
	  for (q in p:ntraits) {
	    if (TraitByTrait[p,q]==1) {
	      Phi.mat[p,q] <- phis[iphi]
	      Phi.mat[q,p] <- Phi.mat[p,q]
	      iphi <- iphi + 1
	    }
	  }
	}
	results <- list(Phi.mat=Phi.mat,phi.log=phi.log)
	return(results)
}
