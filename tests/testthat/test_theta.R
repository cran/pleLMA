context("test_theta.estimates")
library(pleLMA)


test_that("Theta estimates for Rasch, m=1", {

  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)

  model <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
  theta.est <- theta.estimates(inData, model)
 
  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})

test_that("Theta estimates for GPCM, m=1", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
  model <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1)
  theta.est <- theta.estimates(inData, model)

  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})

test_that("Theta estimates for Nominal, m=1", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
  model <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1)
  theta.est <- theta.estimates(inData, model)

  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})


test_that("Theta estimates for Rasch, m=3", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
  inItemTraitAdj <- matrix(c(1,0,0,  1,0,0, 1,0,0,
                             0,1,0,  0,1,0, 0,1,0,
                             0,0,1,  0,0,1, 0,0,1),
                           nrow=9, ncol=3, byrow=TRUE)
  model <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)
  theta.est <- theta.estimates(inData, model)

  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})


test_that("Theta estimates for GPCM, m=3", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
  inItemTraitAdj <- matrix(c(1,0,0,  1,0,0, 1,0,0,
                             0,1,0,  0,1,0, 0,1,0,
                             0,0,1,  0,0,1, 0,0,1),
                           nrow=9, ncol=3, byrow=TRUE)
  model <- ple.lma(inData, model.type="gpcm", inItemTraitAdj, inTraitAdj, tol=1)
  theta.est <- theta.estimates(inData, model)

  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})


test_that("Theta estimates for Nominal, m=3", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
  inItemTraitAdj <- matrix(c(1,0,0,  1,0,0, 1,0,0,
                             0,1,0,  0,1,0, 0,1,0,
                             0,0,1,  0,0,1, 0,0,1),
                           nrow=9, ncol=3, byrow=TRUE)
  model <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj, tol=1)
  theta.est <- theta.estimates(inData, model)

  expect_equal(nrow(theta.est), model$npersons)
  expect_equal(ncol(theta.est), model$ntraits)
})
