# Test helper functions

context("Trimming data frame")



test_that("Trimming results in correct size", {
  data(InstEval)
  trimDat <- trimModelFrame(InstEval)
  expect_more_than(nrow(InstEval), nrow(trimModelFrame(InstEval)))
  expect_equal(nrow(trimDat), 4065)
})
