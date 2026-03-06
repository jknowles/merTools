# Prediction intervals work with new factor levels added part 2

    Code
      predictInterval(glmer3LevSlope, newdata = zNew)
    Condition
      Warning:
      For binomial GLMMs, include.resid.var = TRUE simulates from the
      conditional binomial distribution (n-trial binomial simulation).
      This is the theoretically correct approach.
      To get predictions without residual variance, set include.resid.var = FALSE.
      Warning:
           The following levels of BROOD from newdata 
       -- 100, 101 -- are not in the model data. 
           Currently, predictions for these values are based only on the 
       fixed coefficients and the observation-level error.
    Output
                 fit       upr        lwr
      1   0.34303618 1.0298874 -0.3437464
      2   0.34303618 1.0298877 -0.3437461
      3   0.03703615 0.6510159 -0.6267491
      4   0.06752585 0.7050656 -0.6119811
      5   0.06752641 0.7050647 -0.6119810
      6   0.06752641 0.7050645 -0.6119812
      7   0.06752598 0.7050635 -0.6119814
      8  -0.14685210 0.5119914 -0.8187461
      9  -0.14685250 0.5119904 -0.8187455
      10 -0.14685258 0.5119909 -0.8187455

