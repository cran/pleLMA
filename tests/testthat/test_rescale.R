context("reScaleItem")
library(pleLMA)


test_that("Two Parameter logistic model reScaling of Nominal", {
  data(vocab)
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=10, ncol=1)

  better.sv <- matrix(c(-0.2211536, 0.2211536,
                        -0.6180256, 0.6180256,
                        -0.2175173, 0.2175173,
                        -0.6079710, 0.6079710,
                        -0.4767408, 0.4767408,
                        -0.3758653, 0.3758653,
                        -0.1618439, 0.1618439,
                        -0.2404240, 0.2404240,
                        -0.2126253, 0.2126253,
                        -0.5493119, 0.5493119),
              nrow=10, ncol=2, byrow=TRUE)

  n2pl <- ple.lma(inData=vocab, model.type="nominal", inItemTraitAdj, inTraitAdj,tol=1e-3, starting.sv=better.sv)
  anchor <- c(1,0,0,0,0,0,0,0,0,0)
  new.nu <- reScaleItem(n2pl, anchor)

  expect_equal(sum(new.nu$sNu[1,]), 0)
  expect_equal(sum(new.nu$sNu[1,]**2), 1)

})

test_that("Rescale work for multi-dimensional nominal", {
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

  # Also weak tolerance
  model <- ple.lma(inData, model.type="nominal", inItemTraitAdj, inTraitAdj,
                   tol=1e-01,starting.sv=scores, starting.phi=phi)
  anchor <- c(1,0,0, 0,1,0, 1,0,0)
  new.nu <- reScaleItem(model, anchor)
  
  expect_equal(sum(new.nu$sNu[1,]), 0)
  expect_equal(sum(new.nu$sNu[5,]), 0)
  expect_equal(sum(new.nu$sNu[7,]), 0)
  expect_equal(sum(new.nu$sNu[1,]**2), 1)
  expect_equal(sum(new.nu$sNu[5,]**2), 1)
  expect_equal(sum(new.nu$sNu[7,]**2), 1)
 
})

