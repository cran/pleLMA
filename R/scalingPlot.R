#' Graphs estimated scale values by integers of the LMA (nominal) model
#'
#'  This function plots the estimated item scale values (i.e, nus) by integers
#'  to see shape of scaling of the categories.A linear regression is overlaid
#'  in the plot to help assess linearity. The dashed red line overlaid in the
#'  plot is the linear regression line of the scale values on integers.
#'
#' @param  model.fit   Output from a nominal model
#'
#' @returns plots of estimated scale values by integers
#'
#' @examples
#'
#' #--- some data, 2 items from depression, anxiety and stress scales
#' #    for 250 cases out of possible 1000
#'  data(dass)
#'  inData <- dass[1:250,c("d1", "d2", "a1", "a2", "s1", "s2")]
#'  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
#'  inItemTraitAdj <- matrix(1, nrow=6, ncol=1)
#'  n1 <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-03)
#'  scalingPlot(n1)
#'
#' @export
scalingPlot <- function(model.fit) {
   item.log <- model.fit[[18]]
   nitems <- model.fit[[10]]
   ncat <- model.fit[[11]]
   ItemNames <- model.fit[[5]]

#--- set par back to what it when exit function
   oldpar <- graphics::par(no.readonly = TRUE)    # code line i
   on.exit(graphics::par(oldpar))            # code line i + 1


# --- figure out min and max for axes so that graphs are comparable
	Max <- 0
	Min <- 0
	nus <- matrix(0,nrow=ncat,ncol=nitems)
	for (item in 1:nitems) {
		parameter.record <- item.log[[item]]
		last.parms <- parameter.record[nrow(parameter.record),]
		nu1 <-  -sum(last.parms[(2+ncat):length(last.parms)] )
		nus[,item] <- c(nu1,last.parms[(2+ncat):length(last.parms)])
	if (max(nus) > Max) { Max <- max(nus) }
	if (min(nus) < Min) { Min <- min(nus) }
	}

#  --- for the horizontal axis
	x <- seq(1:ncat)

#  --- and the plots
	graphics::par(mfrow=c(1,2))
	for (item in 1:nitems) {

	ntitle <- paste("Nu: ",ItemNames[item])
	nylab <- expression(nu)

	plot(x,nus[,item],type="b",pch=19,
      main=ntitle,
      ylim=c(Min,Max),
	    ylab= nylab,
	    xlab="Integers")
	    graphics::abline(stats::lm(nus[,item]~x),lty=5,lwd=.7,col="red")
	}
}
