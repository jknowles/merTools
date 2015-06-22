#' @title Extracts random effects
#' @name REextract
#' @description Extracts random effect terms from an lme4 model
#' @param mod a merMod object from the lme4 package
#' @importFrom plyr adply rbind.fill
#' @return a data frame with random effects from lme4 by group id
#' @export
REextract <- function(mod){
  stopifnot(class(mod) == "lmerMod" | class(mod) == "glmerMod")
  out <- lme4::ranef(mod, condVar = TRUE)
  reDims <- length(out)
  tmp.out <- vector("list", reDims)
  for(i in c(1:reDims)){
    tmp.out[[i]] <- out[[i]]
    tmp.out[[i]]$level <- paste0("Level ", i)
    tmp.out[[i]]$level <- as.character(tmp.out[[i]]$level)
    tmp.out[[i]]$groupID <- row.names(out[[i]])
    if(ncol(out[[i]]) > 1){
      tmp.out.se <- plyr::adply(attr(out[[i]], which = "postVar"), c(3),
                              function(x) sqrt(diag(x)))
      colnames(tmp.out.se)[-1] <- paste0(names(out[[i]]), "_se")
      tmp.out.se$X1 <- NULL
      tmp.out[[i]] <- cbind(tmp.out[[i]], tmp.out.se)
    } else {
      tmp.out.se <- sapply(attr(out[[i]], which = "postVar"), sqrt)
      names(tmp.out.se) <- paste0(names(out[[i]]), "_se")
      tmp.out[[i]] <- cbind(tmp.out[[i]], tmp.out.se)
      names(tmp.out[[i]])[4] <-  paste0(names(out[[i]]), "_se")
    }
  }
  dat <- do.call(plyr::rbind.fill, tmp.out)
  # reorg output
  dat[, c("level", "groupID",
          names(dat)[!names(dat) %in% c("level", "groupID")])]
  return(dat)
}

#' @title Simulate random effects from merMod
#' @name REsim
#' @description Simulate random effects from merMod object posterior distributions
#' @param mod a merMod object from the lme4 package
#' @param nsims number of simulations to use
#' @param OR logical, should parameters be converted to odds ratios?
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with distribution of random effect parameters
#' @details Use the Gelman sim technique to build empirical Bayes estimates.
#'  Uses the sim function in the arm package
#' @export
REsim <- function(mod, nsims = 100, OR = FALSE){
  stopifnot(class(mod) == "lmerMod" | class(mod) == "glmerMod")
  mysim <- arm::sim(mod, n.sims = nsims)
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
#' @param mod a merMod object from the lme4 package
#' @param nsims number of simulations to use
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with distribution of fixed effect parameters
#' @details Use the Gelman sim technique to build fixed effect estimates and
#' confidence intervals. Uses the sim function in the arm package
#' @export
FEsim <- function(mod, nsims = 100){
  stopifnot(class(mod) == "lmerMod" | class(mod) == "glmerMod")
  mysim <- arm::sim(mod, n.sims = nsims)
  means <- apply(mysim@fixef, MARGIN = 2, mean)
  medians <- apply(mysim@fixef, MARGIN = 2, median)
  sds <- apply(mysim@fixef, MARGIN =2, sd)
  dat <- data.frame(variable = names(means), mean = means, median = medians,
                    sd = sds, row.names=NULL)
  return(dat)
}

#' @title Estimate the Root Mean Squared Error (RMSE) for a lmerMod
#' @name RMSE.merMod
#' @description Extract the Root Mean Squared Error for a lmerMod object
#' @param mod a lmerMod object from the lme4 package
#' @param scale logical, should the result be returned on the scale of
#' response variable standard deviations?
#' @import lme4
#' @return a numeric which represents the RMSE
#' @export
RMSE.merMod <- function(mod, scale = FALSE){
  stopifnot(class(mod) == "lmerMod")
  # Express RMSE as percentage of dependent variable standard deviation
  dvSD <- sd(mod@frame[, 1])
  RMSE <- sqrt(mean(residuals(mod)^2))
  if(scale == TRUE){
    return(RMSE/dvSD)
  } else{
    return(RMSE)
  }
}

# bootRMSE <- function(data, indx, model){
#   moddat <- data[indx, ]
#   modb <- eval(model@call)
#   RMSE.merMod(modb)
# }
