context("ple_lma_nominal")
library(pleLMA)


test_that("output from fitting uni-dimensional nominal", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)
 # good starting values to speed up test
  scores <- matrix(c(-0.6244714, -0.07818970,  0.17098361, 0.5316775,
                     -0.6136358, -0.12937559,  0.20606888, 0.5369425,
                     -0.6412963, -0.06718668,  0.18456161, 0.5239214,
                     -0.2799779, -0.06605803, -0.01894353, 0.3649794,
                     -0.6211806, -0.21095776, -0.01795778, 0.8500962,
                     -0.3996865, -0.05546513,  0.10737125, 0.3477804,
                     -0.8423171, -0.09575662,  0.31423719, 0.6238365,
                     -0.7782509, -0.13235774,  0.30384290, 0.6067658,
                     -0.9037258, -0.28199819,  0.41075799, 0.7749660),
                   nrow=9, ncol=4, byrow=TRUE)

    model <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj,
                  starting.sv=scores)

  expect_equal(model$npersons, 250)
  expect_equal(model$nitems, 9)
  expect_equal(model$ncat, 4)
  expect_equal(model$nless, 3)
  expect_equal(model$ntraits, 1)
  expect_equal(nrow(model$estimates), 9)
  expect_equal(ncol(model$estimates), 9)
  expect_equal(ncol(model$TraitByTrait), 1)
  expect_equal(nrow(model$TraitByTrait), 1)
  expect_equal(ncol(model$ItemByTrait), 1)
  expect_equal(nrow(model$ItemByTrait), 9)

  expect_is(model, "list")
  expect_is(model$phi.mlogit, "NULL")
  expect_is(model$item.mlogit, "list")
  expect_is(model$Phi.mat, "matrix")
  expect_is(model$mlpl.item, "numeric")
  expect_is(model$mlp.phi, "NULL")

  shows_message("No errors detected in the input")
})


test_that("output from fitting multi-dimensional nominal", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=3, ncol=3)
  inItemTraitAdj <- matrix(c(1,0,0,  1,0,0, 1,0,0,
                             0,1,0,  0,1,0, 0,1,0,
                             0,0,1,  0,0,1, 0,0,1),
                           nrow=9, ncol=3, byrow=TRUE)
  # Good starting values to speed up test
  scores <- matrix(c(-0.9638687, -0.16951993, 0.26351391, 0.8698747,
                 -0.8740057, -0.15198583, 0.26276313, 0.7632284,
                 -1.0227936, -0.20583828, 0.31073510, 0.9178968,
                 -0.4550249, -0.09242915, 0.09653875, 0.4509153,
                 -1.5962158, -0.30439188, 0.16936253, 1.7312451,
                 -0.8293776,  0.01333289, 0.20401866, 0.6120260,
                 -1.1446883, -0.12920615, 0.41535033, 0.8585442,
                 -1.1564133, -0.19243555, 0.41754242, 0.9313064,
                 -1.0652742, -0.44007888, 0.50362027, 1.0017328),
               nrow=9, ncol=4, byrow=TRUE)
    phi <- matrix(c(1.0000000, 0.1754632, 0.2977876,
                    0.1754632, 1.0000000, 0.1934630,
                    0.2977876, 0.1934630, 1.0000000),
                 nrow=3, ncol=3, byrow=TRUE)
   # Also weaker tolerance
   model <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj,
                    tol=1e-01,starting.sv=scores, starting.phi=phi)

  expect_equal(model$npersons, 250)
  expect_equal(model$nitems, 9)
  expect_equal(model$ncat, 4)
  expect_equal(model$nless, 3)
  expect_equal(model$ntraits, 3)
  expect_equal(nrow(model$estimates), 9)
  expect_equal(ncol(model$estimates), 9)
  expect_equal(ncol(model$TraitByTrait), 3)
  expect_equal(nrow(model$TraitByTrait), 3)
  expect_equal(ncol(model$ItemByTrait), 3)
  expect_equal(nrow(model$ItemByTrait), 9)
  expect_equal(sum(diag(model$Phi.mat)),3)

  expect_is(model, "list")
  expect_is(model$item.mlogit, "list")
  expect_is(model$Phi.mat, "matrix")
  expect_is(model$mlpl.item, "numeric")
  expect_is(model$mlpl.phi, "logLik")

  shows_message("No errors detected in the input")
})







