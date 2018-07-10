context("skewness.R")

test_that("skewness works", {
  expect_equal(skewness(septic_patients$age),
               -1.637164,
               tolerance = 0.00001)
  expect_equal(unname(skewness(data.frame(septic_patients$age))),
               -1.637164,
               tolerance = 0.00001)
  expect_equal(skewness(matrix(septic_patients$age)),
               -1.637164,
               tolerance = 0.00001)
})