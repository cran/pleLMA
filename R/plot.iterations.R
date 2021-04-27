#' Plots estimated parameters by iteration for the gpcm and nominal models
#'
#' This is a utility function that plots the estimated item parameters
#' by iterations.  The plots can be used to determine how many iterations
#' are required to get close to final values.  This functions can be used
#' uni- or multi-dimensional gpcm and models.  The number of pages equals
#' the number of items and each page has the plots of marginal effects
#' (left side) and category scale values or alph parameters (right).
#'
#' @param  model.fit    Object from fitting nominal or gpcm model to data
#'
#' @examples
#' 
#' \donttest{
#' data(dass)
#' inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
#' #---   input for uni-dimensional
#' inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#' inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
#'
#' #--- generalized partial credit model
#' g1 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)
#' iterationPlot(g1)
#'
#' #--- nominal response model
#' n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj,inTraitAdj)
#' iterationPlot(n1)
#'
#' 
#' #--- Multidimensional models
#'  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
#'
#' dpress <- matrix(c(1,0,0), nrow=3, ncol=3, byrow=TRUE)
#' anxiety <- matrix(c(0,1,0), nrow=3, ncol=3, byrow=TRUE)
#' stress <- matrix(c(0,0,1), nrow=3, ncol=3, byrow=TRUE)
#' das <- list(dpress, anxiety, stress)
#' inItemTraitAdj <- rbind(das[[1]], das[[2]], das[[3]])
#'
#--- 3 dimensional gpcm
#' g3 <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)
#' iterationPlot(g3)
#'
#--- 3 dimensional nominal
#' n3 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj)
#' iterationPlot(n3)
#' }
#' 
#' @return  Plots of estimated parameters by iteration
#'
#' @export
iterationPlot<- function(model.fit){
  UseMethod("iterationPlot")
}

#' @export

#### for the lambda and nus (or a) over iterations
iterationPlot.list <- function(model.fit) {
	history <- model.fit$item.log
    nitems <- model.fit$nitems
	ncat <- model.fit$ncat
	nless <- model.fit$nless
	ItemNames <- model.fit$ItemNames
	Maxnphi <- model.fit$Maxnphi
	PhiNames <- model.fit$PhiNames
	model.type <- model.fit$model.type
	phi.log <- model.fit$phi.log

	ok.models <- c("nominal", "gpcm")
    "%notin%"<- Negate("%in%")
	if (model.type %notin% ok.models) {
	   stop("Iteration plot only applies to nominal and gpcm  models")
	}

  #--- set par back to what it when exit function
  oldpar <- graphics::par(no.readonly = TRUE)    # code line i
  on.exit(graphics::par(oldpar))            # code line i + 1

 #--- par set in function
 graphics::par(mfrow=c(1,2))


  #--- names of lambda and nu in history
  get.names <- as.data.frame(history[[1]])
  lambda.names <-  names(get.names[3:(2+nless)])
  if (model.type=="nominal") {
     nu.names <-  names(get.names[(2+ncat):(ncat+ncat)])
  } else {
      nu.names <-  names(get.names[(1+ncat)])
  }

  #--- fit overall min and max for y axis values
    histry <- get.names[, lambda.names]
	lam.min<-  min(histry)
    lam.max <- max(histry)
    histry <- get.names[, nu.names]
	  nu.min<-  min(histry)
    nu.max <- max(histry)

  for (item in 2:nitems) {
    histry <- as.data.frame(history[[item]] )
    histry <- histry[2:nrow(histry),]

    lmin <- min(histry[ , lambda.names])
   	if (lmin < lam.min){   lam.min <- lmin }
	lmax <- max(as.matrix(histry[ , lambda.names]))
	if (lmax > lam.max){ lam.max <- lmax}
	nmin <- min(as.matrix(histry[ , nu.names]))
	if (nmin < nu.min) { nu.min <- nmin }
    nmax <- max(as.matrix(histry[ , nu.names]))
	if (nmax > nu.max){ nu.max <- nmax }
	}

  # --- ple item iterations
  for (item in 1:nitems) {

    #--- iteration history for item
    histry <- as.data.frame(history[[item]]                   )
    histry <- histry[2:nrow(histry),]
    iter <- seq(1:nrow(histry))


# Lambda plot
 	# Empty plot to set up
    plot(iter, histry[,1] , type='n',
	    main= paste("Lambda: ","item= ", ItemNames[item]),
		ylab= "Estimated Lambda",
		xlab= "Iteration",
		ylim=c(lam.min,lam.max),
		)

    # Add lines to plot
	for (j in 1:nless) {
		graphics::lines(iter, histry[ , lambda.names[j]], lty=j, lwd=2, col=j)
	  }

# Nu plot
    # Empty plot to set up
	if (model.type == "nominal") {
      plot(iter, histry[,1] , type='n',
	    main= paste("Nu: ","item= ", ItemNames[item]),
		ylab= "Estimated Nu",
		xlab= "Iteration",
		ylim=c(nu.min, nu.max)
		)

	  # Add lines to plot
	    for (j in 1:nless) {
	    graphics::lines(iter, histry[ , nu.names[j]], lty=j, lwd=2, col=j)
	    }
    } else {
      plot(iter, histry[ , nu.names[1]], type='l',  lty=1, lwd=2, col=1,
	    main= paste("a: ","item= ", ItemNames[item]),
		ylab= "Estimated Nu",
		xlab= "Iteration",
		ylim=c(nu.min, nu.max)
		)
	}
  # Next item
  }


# if the number of phis is greater than one...
if (Maxnphi > 1)  {

  #--- names of lambda and nu in history

  phi.log <- as.data.frame(phi.log)
  iter <- 1:nrow(phi.log)
  #--- plot phi.iterations
    for (p in 1:Maxnphi) {
       plot(iter[2:length(iter)], phi.log[2:length(iter) , PhiNames[p]] ,
	     type="l", lwd=2,
         main= PhiNames[p],
         ylim=c(0,1),
         xlab="Iteration",
         ylab=expression(paste("Estimated ", phi)))
  }
}


# end function
}
