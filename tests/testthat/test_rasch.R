context("ple_lma_rasch")
library(pleLMA)


test_that("output from fitting unidimensional Rasch", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
  model <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj)

  expect_equal(model$npersons, 250)
  expect_equal(model$nitems, 9)
  expect_equal(model$ncat, 4)
  expect_equal(model$nless, 3)
  expect_equal(model$ntraits, 1)
  expect_equal(nrow(model$estimates), 9)
  expect_equal(ncol(model$estimates), 8)
  expect_equal(ncol(model$TraitByTrait), 1)
  expect_equal(nrow(model$TraitByTrait), 1)
  expect_equal(ncol(model$ItemByTrait), 1)
  expect_equal(nrow(model$ItemByTrait), 9)

  expect_is(model, "list")
  expect_is(model$phi.mnlogit, "mnlogit")
  expect_is(model$mlpl.phi, "logLik")
  expect_is(model$AIC, "numeric")
  expect_is(model$BIC, "numeric")

  shows_message("No errors detected in the input")

})

test_that("output from fitting multi-dimensional Rasch", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
  inItemTraitAdj <- matrix(c(1,0,0,  1,0,0, 1,0,0,
                             0,1,0,  0,1,0, 0,1,0,
                             0,0,1,  0,0,1, 0,0,1),
                     nrow=9, ncol=3, byrow=TRUE)

  scores <- matrix(c(0,1,2,3),nrow=9,ncol=4, byrow=TRUE)

  model <- ple.lma(inData, model.type="rasch", inItemTraitAdj, inTraitAdj, starting.sv=scores)

  expect_equal(sum(model$estimates[,5:8]-scores), 0)

  expect_equal(model$npersons, 250)
  expect_equal(model$nitems, 9)
  expect_equal(model$ncat, 4)
  expect_equal(model$nless, 3)
  expect_equal(model$ntraits, 3)
  expect_equal(nrow(model$estimates), 9)
  expect_equal(ncol(model$estimates), 8)
  expect_equal(ncol(model$TraitByTrait), 3)
  expect_equal(nrow(model$TraitByTrait), 3)
  expect_equal(ncol(model$ItemByTrait), 3)
  expect_equal(nrow(model$ItemByTrait), 9)

  expect_is(model, "list")
  expect_is(model$phi.mnlogit, "mnlogit")
  expect_is(model$mlpl.phi, "logLik")
  expect_is(model$AIC, "numeric")
  expect_is(model$BIC, "numeric")

  shows_message("No errors detected in the input")

})
