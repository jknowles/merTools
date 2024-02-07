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
      1   0.319305815 1.818685 -1.160422
      2   0.287954733 1.716171 -1.163865
      3   0.054934676 1.448282 -1.360107
      4   0.001274279 1.548561 -1.419932
      5   0.087703116 1.448348 -1.450653
      6   0.034353299 1.563228 -1.428815
      7   0.008824320 1.433124 -1.449002
      8  -0.218615435 1.267652 -1.717512
      9  -0.181019998 1.316339 -1.692194
      10 -0.227707632 1.229211 -1.579833

