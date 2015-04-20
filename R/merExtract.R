#' @title Extracts random effects
#' @name REextract
#' @description Extracts random effect terms from an lme4 model
#' @param mod a merMod object from the lme4 package
#' @importFrom plyr adply
#' @return a data frame with random effects from lme4
#' @export
REextract <- function(mod){
  out <- lme4::ranef(mod, condVar = TRUE)
  out.se <- plyr::adply(attr(out[[1]], which = "postVar"), c(3),
                        function(x) sqrt(diag(x)))
  out.pt <- out[[1]]
  names(out.se)[-1] <- paste0(names(out.pt), "_se")
  newout <- cbind(out.pt, out.se[, -1])
  return(newout)
}

#' @title Simulate random effects from merMod
#' @name REsim
#' @description Simulate random effects from merMod object posterior distributions
#' @param x a merMod object from the lme4 package
#' @param nsims number of simulations to use
#' @param OR logical, should parameters be converted to odds ratios?
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with distribution of random effect parameters
#' @details Use the Gelman sim technique to build empirical Bayes estimates.
#'  Uses the sim function in the arm package
#' @export
REsim <- function(x, nsims, OR = FALSE){
  mysim <- arm::sim(x, n.sims = nsims)
  reDims <- length(mysim@ranef)
  tmp.out <- vector("list", reDims)
  for(i in c(1:reDims)){
    tmp.out[[i]] <- plyr::adply(mysim@ranef[[i]], c(2, 3), plyr::each(c(mean, median, sd)))
    tmp.out[[i]]$level <- paste0("Level ", i)
    tmp.out[[i]]$level <- as.character(tmp.out[[i]]$level)
    tmp.out[[i]]$X1 <- as.character(tmp.out[[i]]$X1)
    tmp.out[[i]]$X2 <- as.character(tmp.out[[i]]$X2)
  }
  dat <- do.call(rbind, tmp.out)
  if(OR == TRUE){
    dat$median <- exp(dat$median)
    dat$mean <- exp(dat$mean)
    dat$sd <- NA # don't know how to do SE of odds ratios currently
    return(dat)
  } else{
    return(dat)
  }
}

#' @title Simulate fixed effects from merMod
#' @name FEsim
#' @description Simulate fixed effects from merMod object posterior distributions
#' @param x a merMod object from the lme4 package
#' @param nsims number of simulations to use
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with distribution of fixed effect parameters
#' @details Use the Gelman sim technique to build fixed effect estimates and
#' confidence intervals. Uses the sim function in the arm package
#' @export
FEsim <- function(x, nsims){
  mysim <- arm::sim(x, n.sims = nsims)
  means <- apply(mysim@fixef, MARGIN = 2, mean)
  medians <- apply(mysim@fixef, MARGIN = 2, median)
  sds <- apply(mysim@fixef, MARGIN =2, sd)
  dat <- data.frame(var = names(means), meanEff = means, medEff = medians,
                    sdEff = sds, row.names=NULL)
  return(dat)
}

# RMSE.merMod <- function(mod){
#   # Express RMSE as percentage of dependent variable standard deviation
#   dvSD <- sd(mod@frame[, 1])
#   RMSE <- sqrt(mean(residuals(mod)^2))
#   return(RMSE/dvSD)
# }
#
# bootRMSE <- function(data, indx, model){
#   moddat <- data[indx, ]
#   modb <- eval(model@call)
#   RMSE.merMod(modb)
# }
