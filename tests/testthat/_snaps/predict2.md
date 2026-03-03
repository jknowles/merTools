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
      1   0.34043520 1.0034767 -0.3276669
      2   0.34043520 1.0034767 -0.3276669
      3   0.03969710 0.6315965 -0.6083550
      4   0.07140516 0.6974340 -0.6097349
      5   0.07140516 0.6974340 -0.6097349
      6   0.07140516 0.6974340 -0.6097349
      7   0.07140516 0.6974340 -0.6097349
      8  -0.14592176 0.4856593 -0.8052807
      9  -0.14592176 0.4856593 -0.8052807
      10 -0.14592176 0.4856593 -0.8052807

