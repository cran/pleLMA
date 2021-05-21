context("y_in_Master")
library(pleLMA)


test_that("y should be factor with 2 levels", {
  data(dass)
  inData <- dass[1:250,c("d1", "d2", "d3", "a1","a2","a3","s1","s2","s3")]
  inTraitAdj  <- matrix(1, nrow=1, ncol=1)
  inItemTraitAdj <- matrix(1, nrow=9, ncol=1)

  i.setup <- set.up(inData, model.type='independence')
  i.master <- i.setup$Master
 
  r.setup <- set.up(inData, model.type='rasch', inTraitAdj,
                  inItemTraitAdj)
  r.master <- r.setup$Master
 
  g.setup <- set.up(inData, model.type='gpcm', inTraitAdj,
                  inItemTraitAdj)
  g.master <- g.setup$Master
 
  n.setup <- set.up(inData, model.type='nominal', inTraitAdj,
                  inItemTraitAdj)
  n.master <- n.setup$Master
 
 expect_equal(dim(table(i.master$y)), 2)
 expect_is(i.master$y, "factor")
 
 expect_equal(dim(table(r.master$y)), 2)
 expect_is(r.master$y, "factor")
 
 expect_equal(dim(table(g.master$y)), 2)
 expect_is(g.master$y, "factor")
 
 expect_equal(dim(table(n.master$y)), 2)
 expect_is(n.master$y, "factor")
 }
 )