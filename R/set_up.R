#' Sets up the data based on input data and model specifications
#'
#' This function sets up the data and sets constants that are essentially
#' the same for all models. This is used within the main wrapper function 
#' `ple.lma', but can also be run independently. If a user wants to run 
#' the functions `fit.independence', `fit.rasch', `fit.gpcm', or `fit.nominal',
#' the set up function should be run prior to using these functions to 
#' create required input. Such an approach can speed up replication studies 
#' because `set.up' would only need to be run once and the response vector 
#' (i.e., named `y') in the Master data frame be replaced by a new one.
#'
#' @param inData       		A person x item Data frame with response patterns
#' @param model.type      Type of model to be fit
#' @param inTraitAdj   		Trait x Trait adjacency matrix (NULL for independence)
#' @param inItemTraitAdj 	Item x Trait adjacency matrix (NULL for independence)
#' @param tol             Tolerence for deteriming convergence (default: 1e-06)
#' @param starting.sv		  Starting category scale values/fixed scores (default: sum equal to zero and sum of squares equal to 1)
#' @param starting.phi    optional: Starting phi matrix (default:  identity matrix)
#'
#' @return PersonByItem   inData (rows are response patterns)
#' @return TraitByTrait   Trait x Trait adjacency matrix
#' @return ItemByTrait    Item x Trait adjacency matrix
#' @return item.by.trait  Need for re-scaling phi.mat
#' @return starting.sv		An item by number of category matrix with starting values for scale values for nominal model and fixed category scores for gpcm and rasch models
#' @return ItemNames      Names of items in inData and PersonByItem
#' @return LambdaName     Short list of lambda names needed for item regressions
#' @return NuName         Short list of nu names names needed for item regressions
#' @return LambdaNames    Long list of lambdas using in Master data set
#' @return NuNames        Long list of nu using in Master data set
#' @return PhiNames       Names of the unique phi parameters
#' @return npersons       Number of individual or persons in data
#' @return nitems         Number of items
#' @return ncat         	Number of categories
#' @return nless          Number of unique lambdas and unique nus
#' @return ntraits        Number of traits
#' @return Maxnphi        Number of phis to estimate
#' @return Nstack         Length of master data set
#' @return pq.mat         An array used to computed (weighted) rest-scores
#' @return Phi.mat        A number of traits x number of traits Phi matrix (defual: the identity matrix)
#' @return Master         Master data set formated for input to to mnlogit
#' @return tol            Tolerence for deteriming convergence
#'
#' @examples
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#'
#'  #--- to set data up for model of independence
#'  ind.setup <- set.up(inData, model.type="independence")
#'
#'  #--- for model specification for uni-dimensional models
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#'  i.setup <- set.up(inData, model.type='independence')
#'  
#'  r.setup <- set.up(inData, model.type='rasch', inTraitAdj,
#'                   inItemTraitAdj)
#'
#'  g.setup <- set.up(inData, model.type='gpcm', inTraitAdj,
#'                   inItemTraitAdj)
#'
#'  n.setup <- set.up(inData, model.type='nominal', inTraitAdj,
#'                   inItemTraitAdj)
#'
#' @export
################################################################################
set.up <- function(inData, model.type, inTraitAdj=NULL, inItemTraitAdj=NULL,
                   tol=NULL, starting.sv=NULL, starting.phi=NULL) {

   #--- CREATE GLOBALS
    # Control algorithm
     if (is.null(tol)) {
       tol <- 1e-6
    } else {
       tol <- tol
    }

  #--- get items names from the data
  ItemNames <- names(inData)

    # Read response patterns into Person by Item data
  PersonByItem <- as.matrix(inData)

     if (!is.null(inTraitAdj)) {
   #--- Read in Trait x Trait Adjaency matrix
       TraitByTrait <- as.matrix(inTraitAdj)
   #--- Read in Item by Trait Adjacency matrix
       ItemByTrait <- as.matrix(inItemTraitAdj)
     }
    # --- Program determined globals based on the data input
    nitems <- length(ItemNames)                  # number of items
    npersons <- nrow(PersonByItem)               # number of persons
    ncat <- dim(table(PersonByItem))             # number of response options
    nless <- ncat-1                              # number of unique lam & nus
    Nstack <- ncat*nitems*npersons               # Number of rows in stacked data


  # --- Names For Master Data
	jcol <- 1
	LambdaNames <- matrix("name",nrow=1,ncol=nitems*nless)
	NuNames <- matrix(0,nrow=1,ncol=nitems)
	for (item in 1:nitems) {
		NuNames[1,item] <- paste("nu",item,sep="")
		for (cat in 1:nless) {
			LambdaNames[1,jcol] <- paste("lam",item,"j",(cat+1), sep="")
			jcol <- jcol + 1
       }
	}
	NuNames <- NuNames
	LambdaNames <- LambdaNames

	# --- Names for ItemFit
	LambdaName <- matrix(,nrow=1,ncol=nless)
	NuName <- matrix(,nrow=1,ncol=nless)
	LambdaNameAll <- matrix(,nrow=1,ncol=nitems*nless)
	for (j in 1:nless) {
	  jj <- j +1
	  LambdaName[1,j] <- paste("lam",jj,sep="")
	  NuName[1,j] <- paste("nu",jj,sep="")
	}
	LambdaName <- LambdaName
	NuName <- NuName

  # Need to create intial MASTER data set
  # ---  Index for Persons
	nrepeats <- ncat*nitems                # Number rows stacked data per person
	personID <- rep(1, nrepeats)            # ID for first person
	for (person in 2:npersons)    {         # ID for rest of them
    personID <- c(personID, rep(person,nrepeats))
	}
	personID <- as.matrix(personID)

  # --- Index for Case
	nrepeats <- ncat;
	CaseIndex <- rep(1,nrepeats)
	for (case in 2:(npersons*nitems)) {
    CaseIndex <- c(CaseIndex, rep(case,nrepeats))
	}
	CaseIndex <- as.matrix(CaseIndex)

  # --- Index for Items:  length ncat*nitem
	ItemIndex <- rep(1,ncat)
	for (item in 2:nitems) {
		ItemIndex <- c(ItemIndex, rep(item,ncat))
	}
	ItemIndex <- as.matrix(ItemIndex)

	ItemID <- ItemIndex
	for (person in 2:npersons) {
		ItemID <- c(ItemID, ItemIndex)
	}

  # --- Index for Categories/response options
	cats <- seq(1:ncat)
	CatIndex <- rep(cats,(nitems*npersons))
	CatIndex <- as.matrix(CatIndex)

  # ---  Block for one item
	block.item <- matrix(0,ncat,nless)     # initial matrix of zeros
	block.item[1,1:nless] <- -1
	for (cat.row in 2:ncat) {
		for (cat.col in 1:nless) {
			if ((cat.row-1) == cat.col)	{
			block.item[cat.row,cat.col] <- 1
			}
		}
	}

	Elambda <- matrix(0,nitems*ncat,nitems*nless)
	irow <- 1
	jcol <- 1
	incat <- ncat
	jnless <- nless
	for (item in 1:nitems) {
		Elambda[irow:incat,jcol:jnless] <- block.item
		irow <- irow + ncat
		jcol <- jcol + nless
		incat <- incat + ncat
	jnless <- jnless + nless
	}
  Elambda.design <- do.call(rbind, replicate(npersons, Elambda, simplify=FALSE))

  # --- Response vector y
	y <- matrix(0,Nstack,1)
	istack <- 1
	for (person in 1:npersons) {
		for (item in 1:nitems) {
			for (cat in 1:ncat) {
				if (cat == PersonByItem[person,item]) {
				y[istack,1] <- 1
				}
		istack <- istack + 1
			}
		}
	}

#--- model specifc parts ---
if (model.type == "independence") {
  # --- Putting it all together for independence
  master <- cbind(personID,CaseIndex,ItemID,CatIndex,y,Elambda.design)
  Master <- as.data.frame(master)
  names(Master) <- c("PersonID", "CaseID", "Item","Category","y",
                     LambdaNames)
} else {

  ntraits <- nrow(TraitByTrait)                # number of latent
  Maxnphi <- 0                                 # number of phi to estimate
  for (i in 1:ntraits) {
    for (j in i:ntraits) {
      if (TraitByTrait[i,j]==1) {
        Maxnphi <- Maxnphi + 1
      }
    }
  }
  Maxnphi <- Maxnphi
  # --- create a (1,nitems) vector with index for latent trait for each item
  item.by.trait <- matrix(0,nrow=1,ncol=nitems)
  for (i in 1:nitems){
    for (j in 1:ntraits) {
      if (ItemByTrait[i,j] == 1) {
        item.by.trait[,i] <- j
      }
    }
  }
  item.by.trait <- item.by.trait  # This is used to impose scaling constraints

   # Phi names
  trait.low <- TraitByTrait & lower.tri(TraitByTrait,diag=T)
  which.low <- which(trait.low, arr.ind=T)
  PhiNames <- paste("phi",apply(which.low, 1, paste, collapse=''), sep='')

  # linear predictor for stacked regressions
  xstack.names <- c(LambdaNames,PhiNames)

  #--- Starting values for the category scale values or setting
  #--- or fixed category scores for rasch and gpcm
  if (is.null(starting.sv) ) {
    starting.sv <- matrix(seq(1:ncat), nrow=nitems, ncol=ncat, byrow=TRUE )
    starting.sv <- starting.sv - mean(starting.sv)         # center
    starting.sv <- starting.sv/sqrt(rowSums(starting.sv**2))    # scale
  } else {
    starting.sv <- starting.sv
  }
  starting.sv <- starting.sv

  # --- Staring values for phis is the identity matrix
	# Control algorithm
	if (is.null(starting.phi)) {
	  Phi.mat <- diag(rep(1,ntraits))
	} else {
	  Phi.mat <- starting.phi
	}


# Now put category scale values in matrix nu
	Nus <- matrix(NA, nrow=npersons*nitems*ncat,ncol=nitems)
  i <- 1
  j <- ncat*nitems
  for (person in 1:npersons) {
		nu.set <- matrix(0,nrow=nitems*ncat,ncol=nitems)
		for (item in 1:nitems) {
			for (category in 1:ncat) {
				if (category==PersonByItem[person,item]) {
				  tmp <- starting.sv[item,category]
				  nu.set[, item] <- tmp
				}
			}
		}
		Nus[i:j,] <- nu.set
		i <- i + ncat*nitems
		j <- j + ncat*nitems
	}


 # --- For summing rest and total scores
	phi.array <- array(0, dim=c(ntraits, ntraits, Maxnphi))
	trait.low <- TraitByTrait & lower.tri(TraitByTrait,diag=T)

 # ---  total number of parameters for latent variables
	which.low <- which(trait.low, arr.ind=T)

 # ---  Will be used to compute rest.scores and totals
	phi.array[cbind(which.low, 1:Maxnphi)] <- 1

 # --- This puts a 1 in the (2,1) element of phi.array[,,2]
	phi.array[cbind(which.low[,2:1], 1:Maxnphi)] <- 1

 # --- pq.mat this is used for rest scores and totals on other traits
	pq.mat <- array(dim=c(nitems, nitems, Maxnphi))
	for (p in 1:Maxnphi){
		pq.mat[,,p] <- ItemByTrait %*% phi.array[,,p] %*% t(ItemByTrait)
		diag(pq.mat[,,p]) <- 0
	}
	pq.mat <- pq.mat

 # --- Putting it all together for rasch, gpcm and nominal models
    master <- cbind(personID,CaseIndex,ItemID,CatIndex,y,Elambda.design,Nus)
    Master <- as.data.frame(master)
    names(Master) <- c("PersonID", "CaseID", "Item","Category","y",
                       LambdaNames,NuNames)
}

# --- these are specifically needed for mnlogit--- may not need
Master$alt <- paste("cat",CatIndex,sep="")
Master$choice <-   with(Master,y == 1)
row.names(Master) <- paste(Master$CaseID,":",Master$alt,sep="")

if (model.type=="independence") {
results <- list(PersonByItem=PersonByItem,
                TraitByTrait= NULL       ,
                ItemByTrait = NULL       ,
                item.by.trait=NULL       ,
                starting.sv = NULL       ,
                ItemNames   =ItemNames   ,
                LambdaName  =LambdaName  ,
                NuName      =NULL        ,
                LambdaNames =LambdaNames ,
                NuNames     =NULL        ,
                PhiNames    =NULL        ,
                npersons    =npersons    ,
                nitems      =nitems      ,
                ncat        =ncat        ,
                nless       =nless       ,
                ntraits     =NULL        ,
                Maxnphi     =NULL        ,
                Nstack      =Nstack      ,
                pq.mat      =NULL        ,
                Phi.mat     =NULL        ,
                Master      =Master
                )
} else if (model.type=="rasch"){
  results <- list(PersonByItem=PersonByItem,
                  TraitByTrait=TraitByTrait,
                  ItemByTrait =ItemByTrait ,
                  item.by.trait=item.by.trait,
                  starting.sv = starting.sv,
                  ItemNames   =ItemNames   ,
                  LambdaName  =LambdaName  ,
                  NuName      =NULL        ,
                  LambdaNames =LambdaNames ,
                  NuNames     =NULL        ,
                  PhiNames    =PhiNames    ,
                  npersons    =npersons    ,
                  nitems      =nitems      ,
                  ncat        =ncat        ,
                  nless       =nless       ,
                  ntraits     =ntraits     ,
                  Maxnphi     =Maxnphi     ,
                  Nstack      =Nstack      ,
                  pq.mat      =pq.mat      ,
                  Phi.mat     =Phi.mat     ,
                  Master      =Master
                )
} else if (model.type=="gpcm") {
  results <- list(PersonByItem=PersonByItem,
                  TraitByTrait=TraitByTrait,
                  ItemByTrait =ItemByTrait ,
                  item.by.trait=item.by.trait,
                  starting.sv = starting.sv,
                  ItemNames   =ItemNames   ,
                  LambdaName  =LambdaName  ,
                  NuName      =NuName      ,
                  LambdaNames =LambdaNames ,
                  NuNames     =NuNames     ,
                  PhiNames    =PhiNames    ,
                  npersons    =npersons    ,
                  nitems      =nitems      ,
                  ncat        =ncat        ,
                  nless       =nless       ,
                  ntraits     =ntraits     ,
                  Maxnphi     =Maxnphi     ,
                  Nstack      =Nstack      ,
                  pq.mat      =pq.mat      ,
                  Phi.mat     =Phi.mat     ,
                  Master      =Master      ,
                  tol         =tol
                  )
} else if (model.type=="nominal"){
  results <- list(PersonByItem=PersonByItem,
                  TraitByTrait=TraitByTrait,
                  ItemByTrait =ItemByTrait ,
                  item.by.trait=item.by.trait,
                  starting.sv = starting.sv,
                  ItemNames   =ItemNames   ,
                  LambdaName  =LambdaName  ,
                  NuName      =NuName      ,
                  LambdaNames =LambdaNames ,
                  NuNames     =NuNames     ,
                  PhiNames    =PhiNames    ,
                  npersons    =npersons    ,
                  nitems      =nitems      ,
                  ncat        =ncat        ,
                  nless       =nless       ,
                  ntraits     =ntraits     ,
                  Maxnphi     =Maxnphi     ,
                  Nstack      =Nstack      ,
                  pq.mat      =pq.mat      ,
                  Phi.mat     =Phi.mat     ,
                  Master      =Master      ,
                  tol         =tol
                  )
}
return <- results
}
