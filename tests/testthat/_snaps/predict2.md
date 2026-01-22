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
      1   0.322506047 1.806669 -1.154082
      2   0.285629074 1.726159 -1.144380
      3   0.055158869 1.443915 -1.348329
      4   0.005072463 1.540368 -1.404610
      5   0.089313444 1.458864 -1.444145
      6   0.030850858 1.554769 -1.402090
      7   0.013920725 1.428768 -1.442980
      8  -0.216355513 1.276630 -1.688061
      9  -0.195913195 1.332680 -1.687023
      10 -0.230126112 1.244005 -1.566616

