#' Fits the nominal model
#'
#' Function estimates the parameters of LMA models where the category scale
#' are estimated. The function can be used to estimate the parameters of the
#' LMA model corresponding the nominal model (for multi-category items) and
#' the 2 parameter logistic model for dichotomous items. The function sets
#' up log object(s) and model formula.  In the case of unidimensional models,
#' the function iterates over item regressions; whereas, for multidimensional
#' models, the function iterates between the item and phi regressions. This
#' function is called from 'ple.lma', but can be run outside of 'ple.lma'.
#'
#' @param  Master		     Master data set in long format
#' @param  Phi.mat		   Matrix of starting values of the association parameters
#' @param  starting.sv   Matrix starting values category scale values
#' @param  pq.mat        Array used compute rest scores and total scores
#' @param  tol           Value used to determine convergence of algorithm
#' @param  PersonByItem  Same as inData (rows are response patterns)
#' @param  TraitByTrait  Same as inTraitAdj (trait x trait adjacency)
#' @param  ItemByTrait   Same as inItemTraitAdj (item x trait adjacency)
#' @param  item.by.trait One dimensional array indicating trait item loads on
#' @param  ItemNames     Names of items in inData (i.e. columns names of
#'                          categorical variables)
#' @param  LambdaNames   Lambda names used in the Master and stacked data frames
#' @param  NuNames       Nu names in Master data frame
#' @param  LambdaName    Lambda names in formula for items
#' @param  NuName        Nu names in formula for item regressions
#' @param  PhiNames      Association parameter names for stacked regression
#' @param  npersons      Number of persons
#' @param  nitems        Number of items
#' @param  ncat          Number of categories per item
#' @param  nless         ncat-1 = number unique lambda and unique nus
#' @param  ntraits       Number of traits
#' @param  Maxnphi       Number of association parametets
#'
#' @return item.log 	 Iteration history of LogLike, lambda, and item parameters
#' @return phi.log		 Iteration history of LogLike, lambdas and phi parameters
#' @return criterion	 Current value of the convergence statistic
#' @return estimates 	 Item x parameter matrix: LogLike, lambda and scale values
#' @return Phi.mat	     Estimated conditional correlation matrix
#' @return fitem		 Formula for item data
#' @return fstack	     Formula for stacked data
#' @return item.mlogit   Summaries from final run of mlogit for item regressions
#' @return phi.mlogit	 Summary from final run of mlogit for stacked regression
#' @return mlpl.item	 Max log pseudo-likelihood function from item regressions
#' @return mlpl.phi      Maximum of log pseudo-likelihood function from stacked regression
#' @return AIC           Akaike information criterion for pseudo-likelihood (smaller is better)
#' @return BIC           Bayesian information criterion for pseudo-likelihood (smaller is better)
#'
#' @examples
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'  #--- unidimensional
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)

#'  s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj,
#'            tol=1e-02)
#'
#'  n1 <- fit.nominal(s$Master, s$Phi.mat, s$starting.sv, s$pq.mat, s$tol,
#'        s$PersonByItem, s$TraitByTrait, s$ItemByTrait, s$item.by.trait,
#'        s$ItemNames, s$LambdaNames,  s$NuNames, s$LambdaName, s$NuName,
#'        s$PhiNames, s$npersons, s$nitems, s$ncat,s$ nless, s$ntraits,
#'        s$Maxnphi)
#'
#' @export
######################################################################################################
fit.nominal <- function(Master, Phi.mat, starting.sv, pq.mat, tol,
                    PersonByItem, TraitByTrait, ItemByTrait, item.by.trait,
                    ItemNames, LambdaNames,  NuNames, LambdaName, NuName,
                    PhiNames, npersons, nitems, ncat,	nless, ntraits, Maxnphi) {

 # --  Object to store results from fitting models to items
	desired_length <- nitems
	item.log <- vector(mode = "list", length = desired_length)
	for (item in 1:nitems) {
		item.log[[item]]<- matrix(0, nrow=1, ncol=(2*nless+2))
	}

  # -- Put starting values first 2 rows (needed if compute in early iterations)
	for (item in 1:nitems) {
		lamholder <- seq(1:nless)
		item.log[[item]] <- c(0 ,0, lamholder, starting.sv[item,2:ncat])
		item.log[[item]] <- rbind(item.log[[item]],
		                          c(0, 0, lamholder, starting.sv[item,2:ncat]) )
	}
  # --- Object to hold iterations from up-dating the phis
    if (ntraits>1) {
     	phiholder <- seq(1:Maxnphi)
	    lamholder <- seq(1:(nitems*nless))
	    phi.log <- matrix(c(0,lamholder,phiholder),nrow=1,
	                      ncol=(1 + nitems*nless + Maxnphi))
	} else {
	    phi.log <- NULL
	}

  # --- Formula for nominal model
  	# up-dating nu
    x.names <- c(LambdaName,NuName)
    fitem <- stats::as.formula(paste("y ~",
                              paste(x.names,collapse="+"),"| 0 | 0",sep=" "))

	# up-dating phi (if needed)
		if (ntraits>1) {
    xstack.names <- c(LambdaNames,PhiNames)
    fstack <- stats::as.formula(paste("y ~",
                         paste(xstack.names,collapse="+"),"| 0 | 0",sep=" "))
		} else {
		    fstack <- NULL
        }


# Fit model
  # --- for uni-dimensional model ---
	if (ntraits==1) {

	criterion <- 10                                   # Any number greater than tol
	while (criterion>tol) {

		ItemLoop.results <- ItemLoop(Master, item.log, Phi.mat=Phi.mat,
			                   PersonByItem=PersonByItem, npersons=npersons,
			                   nitems=nitems, ntraits=ntraits, ncat=ncat,
			                   nless=nless, TraitByTrait=TraitByTrait,
							   pq.mat=pq.mat, LambdaName=LambdaName, NuName=NuName,
							   fitem=fitem)
		Master <- ItemLoop.results$Master
		item.log <- ItemLoop.results$item.log

	  	converge <- convergence.stats(item.log=item.log, nitems=nitems, nless=nless,
		                              LambdaName=LambdaName, NuName=NuName)
			criterion <- abs(converge$criterion.loglike)
			if (criterion > tol) {
			print(paste0(criterion , " > ", tol))
			} else {
				print(paste0("Alogithm has converged: ", criterion," < ",tol))
			}
	}
	}

# --- multidimensional model ---
  if (ntraits>1) {

	criterion <- 10                         # Any number greater than tol

	ItemLoop.results <- ItemLoop(Master, item.log, Phi.mat=Phi.mat,
	                    PersonByItem=PersonByItem, npersons=npersons,
			            nitems=nitems, ntraits=ntraits, ncat=ncat,
	                    nless=nless, TraitByTrait=TraitByTrait,
						pq.mat=pq.mat, LambdaName=LambdaName, NuName=NuName,
						fitem=fitem)
	Master <- ItemLoop.results$Master
	item.log <- ItemLoop.results$item.log

	while (criterion>tol) {

			stack.results <- FitStack(Master, item.log, phi.log, fstack, TraitByTrait,
			                          pq.mat, npersons, nitems, ncat, nless, ntraits,
			                          Maxnphi, PhiNames, LambdaNames)
			Phi.mat <- stack.results[[1]]
			phi.log <- stack.results[[2]]

			scale.results <-Scale(Master, item.log, Phi.mat, PersonByItem,
			                      npersons, nitems, ncat, nless, ntraits,
			                      item.by.trait)
			Master <- scale.results[[1]]
			Phi.mat <- scale.results[[2]]

			ItemLoop.results <- ItemLoop(Master, item.log, Phi.mat=Phi.mat,
			                             PersonByItem=PersonByItem, npersons=npersons,
			                             nitems=nitems, ntraits=ntraits, ncat=ncat,
			                             nless=nless, TraitByTrait=TraitByTrait,
								         pq.mat=pq.mat, LambdaName=LambdaName,
								         NuName=NuName, fitem=fitem)
		    Master <- ItemLoop.results$Master
			  item.log <- ItemLoop.results$item.log

	  converge <- convergence.stats(item.log=item.log, nitems=nitems, nless=nless,
	                                LambdaName=LambdaName, NuName=NuName)
			criterion <- abs(converge$criterion.loglike)
			if (criterion > tol) {
		  print(paste0(criterion , " > ", tol))
			} else {
			print(paste0("Alogithm has converged: ", criterion," < ",tol))
			}
		}
}
	# --- one last set of models to save output from mlogit ---

	desired_length <- nitems
    item.mlogit <- vector(mode = "list", length = desired_length)
    for (item in 1:nitems) {
		DataForItem  <-  ItemData(Master, ItemID=item, Phi.mat=Phi.mat, npersons,
		                          nitems, ntraits, ncat, nless, TraitByTrait,
		                          pq.mat, LambdaName, NuName)
        item.data   <- dfidx::dfidx(DataForItem, choice="y", idx=c("CaseID","alt"))
        item.mlogit[[item]]    <- mlogit::mlogit(fitem, item.data)
 	}

	if (ntraits > 1) {
		   stacked.data <- StackData(Master, item.log, phi.log, pq.mat, npersons,
		                             nitems, ncat, nless, ntraits, Maxnphi,
		                             PhiNames, LambdaNames)
        stacked    <-  dfidx::dfidx(stacked.data, choice="y", idx=c("CaseID","alt"))
        phi.mlogit <- mlogit::mlogit(fstack, stacked)
    } else {
	    phi.mlogit <- NULL
	}

  # --- save LogLike and parmeters to file
	i1 <- item.log[[1]]
	item.save <- i1[nrow(i1),]
	for (item in 2:nitems) {
		i <- item.log[[item]]
		i <- i[nrow(i),]
		item.save <- rbind(item.save,i)
	}
	rownames(item.save) <- ItemNames
	if (ncat > 2) {
	  lam1 <- -rowSums(item.save[,3:(1+ncat)])
	  nu1 <- -rowSums(item.save[,(2+ncat):ncol(item.save)])
	  estimates <- cbind(item.save[,2],lam1,item.save[,3:(1+ncat)],
	                     nu1,item.save[,(2+ncat):ncol(item.save)])
	} else {
     lam1 <- -item.save[,3]
     nu1  <- -item.save[,4]
     estimates <- cbind(item.save[,2],lam1,item.save[,3],nu1,item.save[,4])
	}
	lambdaName<- matrix(NA, nrow=1, ncol=ncat)
	nuName <- matrix(NA,nrow=1,ncol=ncat)
	for (cat in 1:ncat){
	  lambdaName[cat] <- paste("lambda", cat, sep="")
	  nuName[cat] <- paste("nu", cat, sep="")
	}
	colnames(estimates) <- c("loglike", lambdaName, nuName)

	# --- Max of ple function to items and phi
	mlpl.item <- sum(estimates[, 1])
	mlpl.phi  <- phi.mlogit$logLik

	#--- Information criteria
    nparm <- 2*nless*nitems + Maxnphi- ntraits
	AIC <- -1*mlpl.item - nparm
	BIC <- -2*mlpl.item - nparm*log(length(unique(Master$PersonID)))


# --- output from fit.nominal
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
				AIC = AIC[1],
				BIC = BIC[1])
return(results)
}
