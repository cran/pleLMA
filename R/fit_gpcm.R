#' Fits LMA model where category scale values equal a_im * x_j
#'
#' Function estimates the parameters of LMA models with fixed category scores
#' multiplied by an item weight parameter. This function can be used to estimate
#' the LMA model corresponding to is a generalized partial credit model for
#' multi-category items and the 2 parameter logistic model for dichotomous
#' items. The function sets up log objects and model formula.  In the case of
#' unidimensional models, the function iterates over item regressions; whereas,
#' for multidimensional models, the function iterates between the item and phi
#' regressions.  This function is called from 'ple.lma', but can be run outside
#' of 'ple.lma'.
#'
#' @param   Master		     	Master data set in long format
#' @param   Phi.mat		      Matrix of starting values of association parameters
#' @param   PersonByItem  	Person by item matrix of responses (same as inData)
#' @param   TraitByTrait    Trait by trait adjacency matrix (same as inTraitAdj)
#' @param   item.by.trait   Item by trait vector indicating trait item load on (same as inItemTraitAdj)
#' @param   tol			        Criterion used to determine convergence
#' @param   npersons	     	Number of persons
#' @param   nitems			    Number of items
#' @param   ncat			      Number of categories per item
#' @param   nless		      	Number of categories minus 1 (i.e., unique lambdas)
#' @param   ntraits	      	Number of latent traits
#' @param   Maxnphi		      Number of phi parameters to be estimated
#' @param   pq.mat		     	Used to compute rest-scores and totals
#' @param	  starting.sv 	  Fixed category scores
#' @param   ItemNames     	Names of items needed label output
#' @param   LambdaName  	  Names of lambdas needed for formula of the item regressions
#' @param   LambdaNames 	  Names of lambdas needed for formula of the stacked regression
#' @param   PhiNames		    Name of phi parameters (Null for uni-dimensional models)
#'
#' @return 	item.log		  History over iterations of the algorithm for items' log likelihood,
#'                         lambda, and a parameter
#' @return 	phi.log			  History over iterations of the algorithm for log likelihood, lambdas
#'                          nd phi parameters
#' @return 	criterion		  Current value of the convergence statistic which is the maximum of items'
#'                          absolute differences between the current and previous value of the log
#'                          likelihood
#' @return 	estimates   	An item by parameter matrix of estimated item parameter where the
#'                          first column are items' log likelihood
#' @return 	Phi.mat	   	  Estimated matrix of association parameters
#' @return 	fitem		      Formula for item data
#' @return 	fstack	      Formula for stacked data
#' @return 	item.mlogit 	Summary from final run of mlogit for item regressions for each item
#' @return 	phi.mlogit	  Summary from final run mlogit for stacked regression
#' @return 	mlpl.item     Value of maximum of log ple function from fitting items (i.e., sum of logLike)
#' @return 	mlpl.phi	    Value of maximum of log ple function from stacked regression to get phi estimates
#' @return  AIC           Akaike information criterion for pseudo-likelihood (smaller is better)
#' @return  BIC           Bayesian information criterion for pseudo-likelihood (smaller is better)
#'
#' @examples
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "a1","a2","s1","s2")]
#'  #--- unidimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=6, ncol=1)
#'
#' # Need to set up data
#' s <- set.up(inData, model.type='gpcm', inTraitAdj, inItemTraitAdj, tol=1e-03)
#'
#' g <- fit.gpcm(s$Master, s$Phi.mat, s$PersonByItem, s$TraitByTrait,
#'               s$item.by.trait, s$tol, s$npersons, s$nitems, s$ncat,
#'               s$nless, s$ntraits, s$Maxnphi, s$pq.mat, s$starting.sv,
#'               s$ItemNames, s$LambdaName, s$LambdaNames, s$PhiNames)
#'
#' @export

fit.gpcm <- function(Master, Phi.mat, PersonByItem, TraitByTrait, item.by.trait,
                   tol, npersons, nitems, ncat, nless, ntraits, Maxnphi, pq.mat,
                   starting.sv, ItemNames,  LambdaName, LambdaNames, PhiNames) {

  # ---- item history for GPCM -----
	desired_length <- nitems
	item.log <- vector(mode = "list", length = desired_length)
	for (item in 1:nitems) {
		item.log[[item]]<- matrix(0,nrow=1,ncol=(2+ncat))
	}

  # Put starting values in first 2 rows
	for (item in 1:nitems) {
		lamholder <- seq(1:nless)
		item.log[[item]] <- c(0,0,lamholder,1)
		item.log[[item]] <- rbind(item.log[[item]],c(0,0,lamholder,1) )
	}

  # --- phi history  ----
    if (ntraits>1) {
     	phiholder <- seq(1:Maxnphi)
	    lamholder <- seq(1:(nitems*nless))
	    phi.log <- matrix(c(0,lamholder,phiholder),nrow=1,ncol=(1 + nitems*nless + Maxnphi))
	} else {
	    phi.log <- NULL
	}

  # --- formula for GPCM ---
	  gpcm.names <- c(LambdaName,"alpha")
	  fitem <- stats::as.formula(paste("y ~",paste(gpcm.names,collapse="+"),"| 0 | 0",sep=" "))

	# --- formula for up-dating phi (if needed)
	if (ntraits>1) {
			xstack.names <- c(LambdaNames,PhiNames)
			fstack <- stats::as.formula(paste("y ~",paste(xstack.names,collapse="+"),"| 0 | 0",sep=" "))
		} else {
		    fstack <- NULL
        }

  # ---- initial item fit
  	    itemLoop.gpcm <- item.gpcm(Master, item.log, Phi.mat, fitem,
  	                             TraitByTrait,
  	                             PersonByItem, npersons, nitems,
  	                             ncat,  nless, ntraits, Maxnphi, pq.mat,
  	                            starting.sv, LambdaName)
  	    Master <- itemLoop.gpcm$Master
 	    item.log <- itemLoop.gpcm$item.log


  # ---- model fitting loops ----

	if (ntraits==1) {
		criterion <- 10                                   # Any number greater than tol
		while (criterion>tol) {

		    itemLoop.gpcm <- item.gpcm(Master, item.log, Phi.mat, fitem,
		                            TraitByTrait,
		                            PersonByItem, npersons, nitems,
		                            ncat,  nless, ntraits, Maxnphi, pq.mat,
		                            starting.sv, LambdaName)
			Master <- itemLoop.gpcm$Master
			item.log <- itemLoop.gpcm$item.log

			converge <- convergenceGPCM(item.log, nitems, ncat, nless, LambdaName)
			criterion <- abs(converge$criterion.loglike)
			if (criterion > tol) {
				print(paste0(criterion , " > ", tol))
			} else {
				print(paste0("The Alogithm has converged: ", criterion," < ",tol))
			}
		}
	} else {                                      # Multidimensional
	  criterion <- 10                            # Any number greater than tol
	  while (criterion>tol) {

			Stack.results <- fitStackGPCM(Master, item.log, phi.log, fstack,
			                              TraitByTrait, starting.sv, npersons,
			                              nitems, ncat, nless, ntraits, Maxnphi,
			                              pq.mat, LambdaNames, PhiNames)
			Phi.mat <- Stack.results$Phi.mat
			phi.log <- Stack.results$phi.log

			Scale.results <- ScaleGPCM(Master, item.log, Phi.mat, PersonByItem,
			                           npersons, nitems, ncat, nless, ntraits,
			                           starting.sv, item.by.trait)
			Master <- Scale.results[[1]]
			Phi.mat <- Scale.results[[2]]

			itemLoop.gpcm <- item.gpcm(Master, item.log, Phi.mat, fitem,
			                           TraitByTrait,
			                           PersonByItem, npersons, nitems,
			                           ncat,  nless, ntraits, Maxnphi, pq.mat,
			                           starting.sv, LambdaName)
			Master <- itemLoop.gpcm$Master
			item.log <- itemLoop.gpcm$item.log

			converge <- convergenceGPCM(item.log, nitems, ncat, nless, LambdaName)
			criterion <- abs(converge$criterion.loglike)
			if (criterion > tol) {
			print(paste0(criterion , " > ", tol))
			} else {
			print(paste0("The Alogithm has converged: ", criterion," < ",tol))
	    }
	  }
  }

	  # --- one last set of models to save output from mlogit ---
	  desired_length <- nitems
	  item.mlogit <- vector(mode = "list", length = desired_length)
	  for (item in 1:nitems) {
	    DataForItem <- ItemGPCM.data(Master, ItemID=item, Phi.mat, TraitByTrait, pq.mat,
	                             starting.sv,  npersons, nitems, ncat, nless,
	                             ntraits, LambdaName)
        data.fit  <- DataForItem$ItemFit
	      gpcmitem  <- dfidx::dfidx(data.fit, choice="y", idx=c("CaseID","alt"))
        item.mlogit <- mlogit::mlogit(fitem, gpcmitem)
 	  }

	  if (ntraits > 1) {
	    StackedData <- StackDataGPCM(Master, item.log, starting.sv, npersons,
	                                 nitems, ncat, nless, Maxnphi, pq.mat,
	                                 LambdaNames, PhiNames)
	    stacked    <-  dfidx::dfidx(StackedData, choice="y", idx=c("CaseID","alt"))
	    model.fit<- mlogit::mlogit(fstack, stacked)
 	    phi.mlogit <- summary(model.fit)
	  } else {
	      phi.mlogit <- NULL
	  }

  # --- save parameter estimates ---
	i1 <- item.log[[1]]
	item.save <- i1[nrow(i1),]
	for (item in 2:nitems) {
		i <- item.log[[item]]
		i <- i[nrow(i), ]
		item.save <- rbind(item.save,i)
	}
  rownames(item.save) <- ItemNames

	if (ncat > 2) {
        lam1 <- -rowSums(item.save[,3:(1+ncat)])
	    item.save <- cbind(item.save[,2],lam1,item.save[,3:(2+ncat)], starting.sv)
	} else {
       lam1 <- -item.save[,3]
       item.save <- cbind(item.save[,2],lam1,item.save[,3], item.save[,4], starting.sv)
	}
  estimates <- item.save
  lambdaName<- matrix(NA, nrow=1, ncol=ncat)
  xNames <- matrix(NA,nrow=1,ncol=ncat)

  for (cat in 1:ncat){
     lambdaName[cat] <- paste("lambda", cat, sep="")
     xNames[cat] <- paste("x", cat, sep="")
  }
 colnames(estimates) <- c("loglike", lambdaName, "a", xNames)

  # --- Max of ple function to items and phi ---
	mlpl.item <- sum(estimates[, 1])
	mlpl.phi  <- as.numeric(phi.mlogit$logLik)

  # --- Information criteria for pseudo-likelihood
    nparm <- nless*nitems + nitems + Maxnphi- ntraits
    AIC <- -1*mlpl.item - nparm
    BIC <- -2*mlpl.item - nparm*log(length(unique(Master$PersonID)))


	# --- output from fit.gpcm ---
results <- list(item.log=item.log,
                phi.log=phi.log,
                criterion=criterion,
                estimates=estimates,
                Phi.mat=Phi.mat,
                fitem=fitem,
                fstack=fstack,
				        item.mlogit=item.mlogit,
			          phi.mlogit=phi.mlogit,
				        mlpl.item=mlpl.item,
			          mlpl.phi=mlpl.phi,
				        AIC = AIC,
				        BIC = BIC
				      )
return(results)
}
