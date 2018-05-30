# Global variables
utils::globalVariables(c(".shinyMerPar", "sig", "sigma", "Lind", "group",
                         "est"))

#' @importFrom methods as is
#' @importFrom stats AIC as.formula formula logLik median model.matrix na.omit
#' pnorm qnorm quantile residuals rgamma rnorm sd vcov weighted.mean delete.response
#' model.frame na.pass reformulate runif terms getCall
#' @importFrom utils packageVersion
zzz <- function(){
  # Nothing

}


#' Title
#'
#' @param object a merMod object
#' @param correlation optional p value
#' @param use.hessian logical
#' @param ... additional arguments to pass through
#'
#' @return a summary of the object
sum.mm <- function(object,
                           correlation = (p <= getOption("lme4.summary.cor.max")),
                           use.hessian = NULL,
                           ...)
{
  if (length(list(...)) > 0) {
    ## FIXME: need testing code
    warning("additional arguments ignored")
  }
  ## se.calc:
  hess.avail <- (!is.null(h <- object@optinfo$derivs$Hessian) &&
                   nrow(h) > length(getME(object,"theta")))
  if (is.null(use.hessian)) use.hessian <- hess.avail
  if (use.hessian && !hess.avail)
    stop("'use.hessian=TRUE' specified, but Hessian is unavailable")
  resp <- object@resp
  devC <- object@devcomp
  dd <- devC$dims
  ## cmp <- devC$cmp
  useSc <- as.logical(dd[["useSc"]])
  sig <- sigma(object)
  ## REML <- isREML(object)

  famL <- famlink(resp = resp)
  p <- length(coefs <- fixef(object))

  vc <- as.matrix(vcov(object, use.hessian = use.hessian))
  stdError <- sqrt(diag(vc))
  coefs <- cbind("Estimate" = coefs,
                 "Std. Error" = stdError)
  if (p > 0) {
    coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level = 0)
    colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
    if (isGLMM(object)) # FIXME: if "t" above, cannot have "z" here
      coefs <- cbind(coefs, "Pr(>|z|)" =
                       2*pnorm(abs(cf3), lower.tail = FALSE))
  }

  llAIC <- llikAIC(object)
  ## FIXME: You can't count on object@re@flist,
  ##	      nor compute VarCorr() unless is(re, "reTrms"):
  varcor <- VarCorr(object)
  # use S3 class for now
  structure(list(methTitle = methTitle(dd),
                 objClass = class(object),
                 devcomp = devC,
                 isLmer = is(resp, "lmerResp"), useScale = useSc,
                 logLik = llAIC[["logLik"]],
                 family = famL$family, link = famL$link,
                 ngrps = ngrps(object),
                 coefficients = coefs, sigma = sig,
                 vcov = vcov(object, correlation = correlation, sigm = sig),
                 varcor = varcor, # and use formatVC(.) for printing.
                 AICtab = llAIC[["AICtab"]], call = object@call,
                 residuals = residuals(object,"pearson",scaled = TRUE),
                 fitMsgs = fetch.merMod.msgs(object),
                 optinfo = object@optinfo
  ), class = "summary.merMod")
}

#' Find link function family
#'
#' @param object a merMod object
#' @param resp the response vector
#'
#' @return the link function and family
famlink <- function(object, resp = object@resp) {
  if(is(resp, "glmResp"))
    resp$family[c("family", "link")]
  else list(family = NULL, link = NULL)
}


##' Extract all warning msgs from a merMod object
##'
##' @param x a merMod object
fetch.merMod.msgs <- function(x) {
  ## currently only those found with 'X' :
  aX <- attributes(x@pp$X)
  wmsgs <- grep("^msg", names(aX))
  if(any(has.msg <- nchar(Xwmsgs <- unlist(aX[wmsgs])) > 0))
    Xwmsgs[has.msg]
  else
    character()
}
