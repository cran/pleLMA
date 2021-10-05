context("ple_lma_Independence")
library(pleLMA)


test_that("output from fitting independence is correct", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  model <- ple.lma(inData, model.type="independence")

  expect_equal(model$npersons, 250)
  expect_equal(model$nitems, 9)
  expect_equal(model$ncat, 4)
  expect_equal(model$nless, 3)
  expect_equal(nrow(model$estimates), 9)
  expect_equal(ncol(model$estimates), 4)

  expect_is(model, "list")
  expect_is(model$phi.mlogit, "mlogit")
  expect_is(model$mlpl.phi, "numeric")

  shows_message("No errors detected in the input")

})
