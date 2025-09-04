# Prediction intervals work with new factor levels added

    Code
      predictInterval(glmer3LevSlope, newdata = zNew)
    Condition
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
                  fit      upr       lwr
      1   0.319307693 1.818687 -1.160420
      2   0.287956853 1.716173 -1.163863
      3   0.054937049 1.448286 -1.360102
      4   0.001277092 1.548562 -1.419931
      5   0.087706377 1.448351 -1.450651
      6   0.034355538 1.563231 -1.428813
      7   0.008824511 1.433126 -1.449001
      8  -0.218614810 1.267654 -1.717510
      9  -0.181019085 1.316341 -1.692193
      10 -0.227707173 1.229211 -1.579832

