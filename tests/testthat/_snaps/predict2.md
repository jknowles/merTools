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
      1   0.278079625 0.8775483 -0.3470695
      2   0.278079915 0.8775481 -0.3470699
      3   0.001180842 0.6410235 -0.5750499
      4   0.029146030 0.6353424 -0.6320839
      5   0.029146827 0.6353420 -0.6320844
      6   0.029146623 0.6353430 -0.6320840
      7   0.029146526 0.6353431 -0.6320841
      8  -0.164055302 0.4703920 -0.8522081
      9  -0.164055492 0.4703918 -0.8522071
      10 -0.164055141 0.4703921 -0.8522078

