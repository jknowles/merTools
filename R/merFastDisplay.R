#' @title Display model fit summary of merMod objects, fast
#' @name fastdisp
#' @description A faster version of the display function from arm which does
#' not refit the model to extract the deviance
#' @param object a merMod object from the lme4 package
#' @import arm
#' @return A printed summary of a merMod object
#' @export
fastdisp <- function (object, ...)
{
  .local <- function (object, digits = 2, detail = FALSE)
  {
    out <- NULL
    out$call <- object@call
    print(out$call)
    fcoef <- fixef(object)
    useScale <- getME(object, "devcomp")$dims["useSc"]
    corF <- vcov(object)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length(fcoef) > 0) {
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower = FALSE)
        coefs <- cbind(coefs, `z value` = out$z.value,
                       `Pr(>|z|)` = out$p.value)
      }
      else {
        out$t.value <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = out$t.value)
      }
      dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      if (detail) {
        pfround(coefs, digits)
      }
      else {
        pfround(coefs[, 1:2], digits)
      }
    }
    out$coef <- coefs[, "coef.est"]
    out$se <- coefs[, "coef.se"]
    cat("\nError terms:\n")
    vc <- arm:::as.matrix.VarCorr(VarCorr(object), useScale = useScale,
                            digits)
    print(vc[, c(1:2, 4:ncol(vc))], quote = FALSE)
    out$ngrps <- lapply(object@flist, function(x) length(levels(x)))
    is_REML <- isREML(object)
    llik <- logLik(object, REML = is_REML)
    out$AIC <- AIC(llik)
    # out$deviance <- deviance(refitML(object))
    out$n <- getME(object, "devcomp")$dims["n"]
    # Dhat <- -2 * (llik)
    # pD <- out$deviance - Dhat
    # out$DIC <- out$deviance + pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", out$n))
    cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "),
              collapse = "; "))
    cat(sprintf("\nAIC = %g", round(out$AIC, 1)))
    # cat(round(out$DIC, 1))
    # cat("\ndeviance =", fround(out$deviance, 1), "\n")
    if (useScale < 0) {
      out$sigma.hat <- .Call("mer_sigma", object, FALSE,
                             PACKAGE = "lme4")
      cat("overdispersion parameter =", fround(out$sigma.hat,
                                               1), "\n")
    }
    return(invisible(out))
  }
  .local(object, ...)
}
