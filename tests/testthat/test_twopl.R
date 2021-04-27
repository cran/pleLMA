context("ple_2pll")
library(pleLMA)


test_that("Two Parameter Logistic Mode as Nomnial and GPCMl", {
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

  n2pl <- ple.lma(inData=vocab, model.type="nominal", inItemTraitAdj, inTraitAdj,
                   tol=1e-1, starting.sv=better.sv)

  g2pl <- ple.lma(inData=vocab, model.type="gpcm", inItemTraitAdj, inTraitAdj,
                  tol=1e-1)

  gchk <-  round((g2pl$estimates[, 4] * g2pl$estimates[, 5:6]), digits=2)
  nchk <-  round(n2pl$estimates[, 4:5], digits=2)
  chk <- sum(gchk-nchk)

  expect_equal(n2pl$mple.item, g2pl$mple.item)
  expect_equal(chk, 0)

})

