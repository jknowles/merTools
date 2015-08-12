#' @title Display model fit summary of merMod objects, fast
#' @name fastdisp
#' @description A faster version of the \code{arm::\link[arm]{display}} function that
#' is quicker because it does not refit the model to extract the deviance
#' @param merMod a merMod object from the lme4 package
#' @param ... additional arguments to pass to \code{arm::\link[arm]{display}}
#' including number of digits
#' @details The time saving is only noticeable for large, time-consuming (g)lmer
#' fits.
#' @import arm
#' @return A printed summary of a merMod object
#' @examples
#' \dontrun{
#' #Compare the time for displaying this modest model
#' require(arm)
#' m1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
#' system.time(display(m1))
#' system.time(fastdisp(m1))
#' }
#' @seealso \code{\link[arm]{display}}
#' @export
fastdisp <- function (merMod, ...)
{
  .local <- function (merMod, digits = 2, detail = FALSE)
  {
    out <- NULL
    out$call <- merMod@call
    print(out$call)
    fcoef <- fixef(merMod)
    useScale <- getME(merMod, "devcomp")$dims["useSc"]
    corF <- vcov(merMod)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length(fcoef) > 0) {
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
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
    vc <- easyVarCorr(VarCorr(merMod), useScale = useScale,
                            digits)
    print(vc[, c(1:2, 4:ncol(vc))], quote = FALSE)
    out$ngrps <- lapply(merMod@flist, function(x) length(levels(x)))
    is_REML <- isREML(merMod)
    llik <- logLik(merMod, REML = is_REML)
    out$AIC <- AIC(llik)
    # out$deviance <- deviance(refitML(merMod))
    out$n <- getME(merMod, "devcomp")$dims["n"]
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
      out$sigma.hat <- lme4::sigma(merMod)
      cat("overdispersion parameter =", fround(out$sigma.hat,
                                               1), "\n")
    }
    return(invisible(out))
  }
  .local(merMod, ...)
}
