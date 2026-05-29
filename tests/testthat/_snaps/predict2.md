# Prediction intervals work with new factor levels added part 2

    Code
      predictInterval(glmer3LevSlope, newdata = zNew)
    Condition
      Warning:
      For binomial GLMMs, include.resid.var = TRUE simulates from the
      conditional binomial distribution (n-trial binomial simulation).
      This is the theoretically correct approach.
      To get predictions without residual variance, set include.resid.var = FALSE.
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning in `chol.default()`:
      the matrix is either rank-deficient or not positive definite
      Warning:
           The following levels of BROOD from newdata 
       -- 100, 101 -- are not in the model data. 
           Currently, predictions for these values are based only on the 
       fixed coefficients and the observation-level error.
    Output
                  fit       upr        lwr
      1   0.286184362 0.8728483 -0.3369354
      2   0.286184362 0.8728483 -0.3369354
      3  -0.004346961 0.6208335 -0.5510659
      4   0.032606238 0.6251944 -0.6043012
      5   0.032606238 0.6251944 -0.6043012
      6   0.032606238 0.6251944 -0.6043012
      7   0.032606238 0.6251944 -0.6043012
      8  -0.161448918 0.4544298 -0.8194496
      9  -0.161448918 0.4544298 -0.8194496
      10 -0.161448918 0.4544298 -0.8194496

