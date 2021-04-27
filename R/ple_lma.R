#' Main function for estimating parameters of LMA models
#'
#' This function is a wrapper function that checks for errors in the input
#'  (i.e., `error.check'),  sets up required objects and data (i.e., `set.up'),
#'  calls function to fit specified model (either `fit.independence',
#'  fit.rasch', `fit.gpcm', or  `fit.nomial'), and outputs an extensive
#'  list of details and results. The required input for all models consist
#'  of a data frame where elements are consecutive integers (1, 2, ...)
#'  indicating the category chosen by each case/individual (rows) for
#'  each variable (columns), and the model type.  For the LMA models that
#'  correspond to item response theory models require an Item x Trait
#'  adjacency matrix (`inItemTraitAdj') and a Trait x Trait adjacency matrix
#'  (`inTraitAdj').  Optional input  include the tolerance  value (`tol')
#'  which is used to for determine whether the pseudo-likelihood algorithm
#'  has converged for a gpcm or nominal model default=1e-6).  Additional
#'  optional input (`starting.sv') is an item by category a matrix of
#'  starting scale values for the nominal model or fixed category scores
#'  for the gpcm and rasch models. The default category scale values/scores
#'  are eqaully spaced, centered at zero, and the sum of squared values
#'  equals 1. The final optional input is a trait x trait (`starting.phi')
#'  matrix of starting values for the association parameter matrix
#'  (default= identity matrix).
#'
#' @param inData          A person x item data matrix or data frame with elements equal to reponse options choosen by an individual.
#' @param inItemTraitAdj  An Item x Trait adjacency matrix indicating what trait an item loads on.
#' @param inTraitAdj      A Trait x Trait adjacency matrix indicating relationship among traits.
#' @param model.type   	  Model to be fit (nominal, gpcm, rasch, independence) to data.
#' @param tol         	  Convergence criterion, default = 1e-6
#' @param starting.sv	    Starting category scale values/fixed scores
#' @param starting.phi    Starting matrix of phi parameters (i.e., conditional covariance matrix)
#'
#' @return model.type     The model (nominal, gpcm, rash, or independence) that was fit to data
#' @return TraitByTrait   The Trait x Trait adjacency matrix used.
#' @return ItemByTrait    The Item x Trait adjacency matrix.
#' @return item.by.trait  One dimensional version of ItemByTrait that gives the number of trait.
#' @return ItemNames      Names of items in inData
#' @return PhiNames       Names of the association parameters (i.e., phi)
#' @return formula.item   Formula used to up-date item parameters via item regressions.
#' @return formula.phi    Formula used to up-date association parameters via stacked regression.
#' @return npersons       Number of persons in data set.
#' @return nitems         Number of items.
#' @return ncat           Number of categories per item.
#' @return nless          Number of unique marginal effects & unique scale values.
#' @return Maxnphi        Number of association parameters estimated.
#' @return ntraits        Number of traits.
#' @return starting.sv    Starting scale values for nominal model or fixed scores for rasch or gpcm.
#' @return tol            Used to determine convergence default= 1e-7
#' @return criterion      Final value criterion at convergence
#' @return item.log       Item iteration history plus maximum of the LogLike for each item
#' @return phi.log        Assocation parameter iteration history
#' @return estimates      Item x Parameter matrix where 1st column is max LogLike for each item and remaining columns are item parameter estimate
#' @return Phi.mat        Estimated conditional correlation matrix
#' @return item.mnlogit   Output from final mnlogit fit to items
#' @return phi.mnlogit    Output form final mnlogit fit to stacked data
#' @return mlpl.item      Max Log(pseudo-likelihood) function from item models (i.e. sum of first column of estimates)
#' @return mlpl.phi       Max Log(pseudo-likelihood) function from stacked regression(s).
#' @return AIC            Akaike information criterion for pseudo-likelihood (smaller is better)
#' @return BIC            Bayesian information criterion for pseudo-likelihood (smaller is better)
#'
#' @examples
#' #--- some data, 3 items from dpression, anxiety and stress scales
#' #    and only 250 cases out of possible 1000
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#' #---  log-linear model of independence
#'  ind <- ple.lma(inData, model.type="independence")
#'
#' #---   input for uni-dimensional
#'   inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'   inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#' #--- rasch family
#'  r1 <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
#'
#' #---  rasch with alternative scores
#'  scores <- matrix(c(0,1,2,3),nrow=9,ncol=4,byrow=TRUE)
#'  r1b <- ple.lma(inData, model.type="rasch", inItemTraitAdj,
#'                 inTraitAdj, starting.sv=scores)
#'
#' \donttest{
#' #--- generalized partial credit model
#' g1 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)
#'
#' #--- gpcm with alternative scores
#' scores <- matrix(c(0,1,2,3),nrow=9,ncol=4,byrow=TRUE)
#' g1b <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj, starting.sv=scores)
#'
#' #--- nominal response model
#' n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj,inTraitAdj)
#'
#' #--- re-run nominal model with input starting values and phi
#' #    and setting stronger convergnce criterion.
#' sv <- n1$estimates[, 6:9]
#' phi <- n1$Phi.mat
#' n1b <- ple.lma(inData, model.type="nominal", inItemTraitAdj,
#'                inTraitAdj, starting.sv=sv, starting.phi=phi, tol=1e-8)
#' }
#' 
#' 
#' #--- Multidimensional models
#' #--- re-define inTraitAdj and inItemTraitAdj for 3 dimensional models
#'   inTraitAdj  <- matrix(1, nrow=3, ncol=3)
#'
#'   dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
#'   anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
#'   stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
#'   das <- list(dpress, anxiety, stress)
#'   inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])
#'
#' #--- 3 dimensional rasch
#'   r3 <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
#'
#' \donttest{
#' #--- 3 dimensional gpcm
#'   g3 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)
#'
#' #--- 3 dimensional nominal
#'   n3 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj)
#' 
#' 
#' #--- 2 parameter logistic IRT model fit to responses to
#' #    10 dichotomous Vocabulary items from from 2018 GSS
#' #    by 1309 respondents
#'   data(vocab)
#'   inItemTraitAdj <- matrix(1, nrow=10, ncol=1)
#'   inTraitAdj <- matrix(1, nrow=1, ncol=1)
#'
#' #--- 2 pl as a gpcm model
#'   g.2pl <- ple.lma(inData=vocab, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1e-03)
#'
#' #--- 2 pl as a nominal model
#'   n.2pl <- ple.lma(inData=vocab, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
#' }
#'
#' @export
###############################################################################
ple.lma <- function(inData, model.type, inItemTraitAdj=NULL, inTraitAdj=NULL,
                    tol=NULL, starting.sv=NULL, starting.phi=NULL) {


  #--- simple check for errors in the data input ---
  error.check(inData, model.type, inTraitAdj, inItemTraitAdj)

  #--- control algorithm stopped criteria
  if (is.null(tol)) {
    tol <- 1e-6
  } else {
    tol <- tol
  }

  #--- set up master data set, constants and other important things that are
  #    passed into functions
  setup <- set.up(inData, model.type, inTraitAdj, inItemTraitAdj, tol,
                  starting.sv, starting.phi)
  tol          <- setup$tol
  PersonByItem <- setup$PersonByItem
  TraitByTrait <- setup$TraitByTrait
  ItemByTrait  <- setup$ItemByTrait
  item.by.trait<- setup$item.by.trait
  starting.sv  <- setup$starting.sv
  ItemNames    <- setup$ItemNames
  LambdaName   <- setup$LambdaName
  NuName       <- setup$NuName
  LambdaNames  <- setup$LambdaNames
  NuNames      <- setup$NuNames
  PhiNames     <- setup$PhiNames
  npersons     <- setup$npersons
  nitems       <- setup$nitems
  ncat         <- setup$ncat
  nless        <- setup$nless
  ntraits      <- setup$ntraits
  Maxnphi      <- setup$Maxnphi
  Nstack       <- setup$Nstack
  pq.mat       <- setup$pq.mat
  Phi.mat      <- setup$Phi.mat
  Master       <- setup$Master

  message("Basic set up is complete")

#####################################################################
# Ready to go to modeling fitting                                   #
#####################################################################
if (model.type == "nominal") {

	    model.results <- fit.nominal(Master, Phi.mat, starting.sv, pq.mat, tol,
	                     PersonByItem, TraitByTrait, ItemByTrait, item.by.trait,
	                     ItemNames, LambdaNames, NuNames, LambdaName, NuName,
	                     PhiNames, npersons, nitems, ncat, nless, ntraits,
	                     Maxnphi)

} else if (model.type== "gpcm") {

	    model.results <- fit.gpcm(Master, Phi.mat, PersonByItem, TraitByTrait,
	                              item.by.trait, tol, npersons, nitems, ncat,
	                              nless, ntraits, Maxnphi, pq.mat, starting.sv,
	                              ItemNames, LambdaName, LambdaNames, PhiNames)

} else if (model.type=="rasch") {

	      model.results  <- fit.rasch(Master, npersons, nitems, ncat, nless,
	                                  Maxnphi, pq.mat,	starting.sv, LambdaNames,
	                                  PhiNames, ItemNames, LambdaName)

 } else  if (model.type == "independence") {
      model.results <- fit.independence(Master, LambdaNames, LambdaName,
                                        ItemNames)
 }

all.results <- list(model.type   = model.type,
              TraitByTrait = TraitByTrait,
			        ItemByTrait  = ItemByTrait,
			        item.by.trait= item.by.trait,
	            ItemNames    = ItemNames,
			        PhiNames     = PhiNames,
		          formula.item = model.results$fitem,
	            formula.phi  = model.results$fstack,
              npersons     = npersons,
	            nitems       = nitems,
		          ncat         = ncat,
					    nless        = nless,
					    Maxnphi      = Maxnphi,
	            ntraits      = ntraits,
		          starting.sv  = starting.sv,
		          tol          = tol,
		          criterion    = model.results$criterion,
		          item.log     = model.results$item.log,
	            phi.log      = model.results$phi.log,
		          estimates    = model.results$estimates,
		          Phi.mat      = model.results$Phi.mat,
		          item.mnlogit = model.results$item.mnlogit,
		          phi.mnlogit  = model.results$phi.mnlogit,
	            mlpl.item    = model.results$mlpl.item,
	            mlpl.phi     = model.results$mlpl.phi,
					    AIC          = model.results$AIC,
					    BIC          = model.results$BIC)

return(all.results)
}
