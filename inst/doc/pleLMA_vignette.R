## -----------------------------------------------------------------------------
library(pleLMA)

## -----------------------------------------------------------------------------
data(dass)

## ---- eval=FALSE--------------------------------------------------------------
#  ?dass

## -----------------------------------------------------------------------------
data(dass)
items.to.use <- c("d1","d2","d3","a1","a2","a3","s1","s2","s3")
inData <- dass[1:250,(items.to.use)]
head(inData)

## -----------------------------------------------------------------------------
#--- Log-linear model of Independence
ind <- ple.lma(inData, model.type="independence")

## -----------------------------------------------------------------------------
(inTraitAdj <- matrix(1, nrow=1 ,ncol=1))

## -----------------------------------------------------------------------------
(inItemTraitAdj <- matrix(1, nrow=ncol(inData), ncol=1) )

## -----------------------------------------------------------------------------
#--- Model in the rasch family
r1 <-  ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)

## -----------------------------------------------------------------------------
#--- Generalized partial credit model
g1 <-  ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)

## -----------------------------------------------------------------------------
#--- Nominal response model
n1 <-  ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj)

## -----------------------------------------------------------------------------
xj <- matrix(c(0, 1, 2, 5), nrow=9, ncol=4, byrow=TRUE)
g1b <- ple.lma(inData, inItemTraitAdj, inTraitAdj, model.type="gpcm", starting.sv=xj)

## -----------------------------------------------------------------------------
n1$estimates
sv <- n1$estimates[, 6:9]
sigma <- n1$Phi.mat
n1.alt <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj,
                  starting.sv = sv,  starting.phi= sigma)

## ---- eval=FALSE--------------------------------------------------------------
#  n1.summary <- lma.summary(n1)

## -----------------------------------------------------------------------------
noquote(lma.summary(n1)$report)

## -----------------------------------------------------------------------------
lma.summary(n1)$TraitByTrait
lma.summary(n1)$ItemByTrait

## -----------------------------------------------------------------------------
#--- item by log likelihoods, lambdas, and nus
lma.summary(n1)$estimates   
#--- sigma_1^2
lma.summary(n1)$phi

## -----------------------------------------------------------------------------
g1$estimates

## ---- eval=FALSE--------------------------------------------------------------
#  n1$item.log[[1]]

## ---- eval=FALSE--------------------------------------------------------------
#  iterationPlot(n1)

## -----------------------------------------------------------------------------
s <- set.up(inData, model.type='nominal', inTraitAdj, inItemTraitAdj)

convergence.stats(n1$item.log, n1$nitems, n1$nless, s$LambdaName, s$NuName)

## ---- eval=FALSE--------------------------------------------------------------
#  s <- set.up(inData, model.type='gpcm', inTraitAdj, inItemTraitAdj)
#  
#  convergenceGPCM(g1$item.log, g1$nitems, g1$ncat, g1$nless, s$LambdaName)

## ---- eval=FALSE--------------------------------------------------------------
#  scalingPlot(n1)

## -----------------------------------------------------------------------------
anchor <- matrix(0, nrow=1, ncol=9)
anchor[1,1] <- 1
    
reScaleItem(n1, anchor=anchor) 

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  theta.r1 <- theta.estimates(inData, r1)
#  theta.g1 <- theta.estimates(inData, g1)
#  theta.n1 <- theta.estimates(inData, n1)

## -----------------------------------------------------------------------------
(inTraitAdj <- matrix(1, nrow=3 ,ncol=3))

## -----------------------------------------------------------------------------
d <- matrix(c(1, 0, 0),nrow=3,ncol=3,byrow=TRUE)
a <- matrix(c(0, 1, 0),nrow=3,ncol=3,byrow=TRUE)
s <- matrix(c(0, 0, 1),nrow=3,ncol=3,byrow=TRUE)
das <- list(d, a, s)
(inItemTraitAdj  <- rbind(das[[1]], das[[2]], das[[3]]))

## ---- echo=TRUE, results='hide'-----------------------------------------------
r3 <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
  
g2 <-  ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj)

n3 <-  ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj)

## -----------------------------------------------------------------------------
noquote(lma.summary(n3)$report) 

## -----------------------------------------------------------------------------
n3$Phi.mat

## -----------------------------------------------------------------------------
# the full data set
inData <- dass

# A (3 x 3) trait by trait adjacency matrix
inTraitAdj <- matrix(c(1,1,1, 1,1,1, 1,1,1), nrow=3 ,ncol=3)

# A (42 x 3) item by trait adjacency matrix
d <- matrix(c(1, 0, 0),nrow=14,ncol=3,byrow=TRUE)
a <- matrix(c(0, 1, 0),nrow=13,ncol=3,byrow=TRUE)
s <- matrix(c(0, 0, 1),nrow=15,ncol=3,byrow=TRUE)
das <- list(d, a, s)
inItemTraitAdj  <- rbind(das[[1]], das[[2]], das[[3]])

## -----------------------------------------------------------------------------
data(vocab)
inItemTraitAdj <- matrix(1, nrow=ncol(vocab), ncol=1)
inTraitAdj <- matrix(1, nrow=1, ncol=1)

## ---- echo=TRUE, results='hide'-----------------------------------------------
#--- 2 pl as a gpcm model
  g.2pl <- ple.lma(inData=vocab, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1e-04)

#--- 2 pl as a nominal model
  n.2pl <- ple.lma(inData=vocab, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1e-04)

## -----------------------------------------------------------------------------
g.2pl$estimates

n.2pl$estimates

## -----------------------------------------------------------------------------
g.2pl$estimates[, 4] * g.2pl$estimates[, 5:6]

