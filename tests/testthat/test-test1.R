context("running linear_regression function")

test_that("running test case #1", {
  #define simple predictor and target
  Y = matrix(c(5.6, 7.9, 10.8), ncol=1)
  X = matrix(c(1, 2, 3), ncol=1)
  #run liner regression
  reg = linear_regression(X, Y)

  #define reference values
  beta0_ref = 2.9
  beta1_ref = 2.6
  beta0 = reg[1,1]
  beta1 = reg[2,1]

  expect_equal(beta0, beta0_ref)
  expect_equal(beta1, beta1_ref)
  }
)
