# Global variables
utils::globalVariables(c(".shinyMerPar", "sig"))

# Namespace fixes

#' @importFrom methods as
#' @importFrom stats AIC as.formula formula logLik median model.matrix na.omit
#' pnorm qnorm quantile residuals rgamma rnorm sd vcov weighted.mean
zzz <- function(){
  # Nothing

}
