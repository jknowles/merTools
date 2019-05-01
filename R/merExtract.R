#' @title Extracts random effects
#' @name REextract
#' @description Extracts random effect terms from an lme4 model
#' @param merMod a merMod object from the lme4 package
#' @return a data frame with the following columns
#' \describe{
#'   \item{groupFctr}{The name of the grouping factor associated with the random effects}
#'   \item{groupID}{The level of the grouping factor associated with the random effects}
#'   \item{'term'}{One column per random effect, the name is derived from the merMod}
#'   \item{'term'_se}{One column per random effect, the name is derived from the merMod}
#' }
#' @examples
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' rfx <- REextract(m2)
#' #Note the column names
#' head(rfx)
#' @export
REextract <- function(merMod){
  stopifnot(inherits(merMod, "merMod"))
  
  out <- lme4::ranef(merMod, condVar = TRUE)
  lvlNames <- names(out)
  reDims <- length(out)
  tmp.out <- vector("list", reDims)
  for(i in c(1:reDims)){
    tmp.out[[i]] <- out[[i]]
    tmp.out[[i]]$groupFctr <- lvlNames[i]
    tmp.out[[i]]$groupID <- row.names(out[[i]])
    if(ncol(out[[i]]) > 1){
      tmp.out.se <- apply(attr(out[[i]], which = "postVar"), 3,
                          function(x) sqrt(diag(x)))
      tmp.out.se <- as.data.frame(t(tmp.out.se))
      colnames(tmp.out.se) <- paste0(names(out[[i]]), "_se")
      tmp.out[[i]] <- cbind(tmp.out[[i]], tmp.out.se)
    } else {
      tmp.out.se <- sapply(attr(out[[i]], which = "postVar"), sqrt)
      names(tmp.out.se) <- paste0(names(out[[i]]), "_se")
      tmp.out[[i]] <- cbind(tmp.out[[i]], tmp.out.se)
      names(tmp.out[[i]])[4] <-  paste0(names(out[[i]]), "_se")
    }
  }
  dat <- dplyr::bind_rows(tmp.out)
  # reorg output
  dat <- dat[, c("groupFctr", "groupID",
          names(dat)[!names(dat) %in% c("groupFctr", "groupID")])]
  return(dat)
}

#' Simulate random effects from merMod
#' \code{REsim} simulates random effects from merMod object posterior distributions
#' @param merMod a merMod object from the lme4 package
#' @param n.sims number of simulations to use
#' @param oddsRatio logical, should parameters be converted to odds ratios?
#' @param seed numeric, optional argument to set seed for simulations
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with the following columns
#' \describe{
#'   \item{\code{groupFctr}}{Name of the grouping factor}
#'   \item{\code{groupID}}{Level of the grouping factor}
#'   \item{\code{term}}{Name of random term (intercept/coefficient)}
#'   \item{\code{mean}}{Mean of the simulations}
#'   \item{\code{median}}{Median of the simulations}
#'   \item{\code{sd}}{Standard deviation of the simulations, \code{NA} if \code{oddsRatio=TRUE}}
#' }
#' @details Use the Gelman sim technique to build empirical Bayes estimates.
#'  Uses the sim function in the arm package
#' @examples
#' require(lme4)
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' re2 <- REsim(m2, 25)
#' head(re2)
#' @export
REsim <- function(merMod, n.sims = 200, oddsRatio = FALSE, seed=NULL){
  stopifnot(inherits(merMod, "merMod"))
  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

  mysim <- arm::sim(merMod, n.sims = n.sims)
  reDims <- length(mysim@ranef)
  tmp.out <- vector("list", reDims)
  names(tmp.out) <- names(mysim@ranef)
  for(i in c(1:reDims)){
    zed <- apply(mysim@ranef[[i]], c(2, 3),
                 function(x) as.data.frame(x) %>% dplyr::summarise_all(.funs = c("mean", "median", "sd")))
    zed <- bind_rows(zed)
    zed$X1 <- rep(dimnames(mysim@ranef[[i]])[[2]], length(dimnames(mysim@ranef[[i]])[[3]]))
    zed$X2 <- rep(dimnames(mysim@ranef[[i]])[[3]], each = length(dimnames(mysim@ranef[[i]])[[2]]))
    tmp.out[[i]] <- zed; rm(zed)
    tmp.out[[i]]$groupFctr <- names(tmp.out)[i]
    tmp.out[[i]]$X1 <- as.character(tmp.out[[i]]$X1)
    tmp.out[[i]]$X2 <- as.character(tmp.out[[i]]$X2)
  }
  dat <- do.call(rbind, tmp.out)
  dat$groupID <- dat$X1; dat$X1 <- NULL
  dat$term <- dat$X2; dat$X2 <- NULL
  dat <- dat[, c("groupFctr", "groupID", "term", "mean", "median", "sd")]
  rownames(dat) <- NULL
  if(oddsRatio == TRUE){
    dat$median <- exp(dat$median)
    dat$mean <- exp(dat$mean)
    dat$sd <- NA # don't know how to do SE of odds ratios currently
    return(dat)
  } else{
    return(dat)
  }
}

#' Simulate fixed effects from merMod
#' \code{FEsim} simulates fixed effects from merMod object posterior distributions
#' @param merMod a merMod object from the lme4 package
#' @param n.sims number of simulations to use
#' @param oddsRatio logical, should parameters be converted to odds ratios?
#' @param seed numeric, optional argument to set seed for simulations
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with the following columns
#' \describe{
#'   \item{\code{term}}{Name of fixed term (intercept/coefficient)}
#'   \item{\code{mean}}{Mean of the simulations}
#'   \item{\code{median}}{Median of the simulations}
#'   \item{\code{sd}}{Standard deviation of the simulations, \code{NA} if \code{oddsRatio=TRUE}}
#' }
#' @details Use the Gelman sim technique to build fixed effect estimates and
#' confidence intervals. Uses the sim function in the arm package
#' @examples
#' require(lme4)
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fe2 <- FEsim(m2, 25)
#' head(fe2)
#' @export
FEsim <- function(merMod, n.sims = 200, oddsRatio=FALSE, seed=NULL){
  stopifnot(inherits(merMod, "merMod"))
  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

  mysim <- arm::sim(merMod, n.sims = n.sims)
  means <- apply(mysim@fixef, MARGIN = 2, mean)
  medians <- apply(mysim@fixef, MARGIN = 2, median)
  sds <- apply(mysim@fixef, MARGIN =2, sd)
  dat <- data.frame(term = names(means), mean = means,
                    median = medians,
                    sd = sds, row.names=NULL)
  if(oddsRatio == TRUE){
    dat$median <- exp(dat$median)
    dat$mean <- exp(dat$mean)
    dat$sd <- NA # don't know how to do SE of odds ratios currently
    return(dat)
  } else{
    return(dat)
  }
}

#' @title Estimate the Root Mean Squared Error (RMSE) for a lmerMod
#' @name RMSE.merMod
#' @description Extract the Root Mean Squared Error for a lmerMod object
#' @param merMod a lmerMod object from the lme4 package
#' @param scale logical, should the result be returned on the scale of
#' response variable standard deviations?
#' @import lme4
#' @return a numeric which represents the RMSE
#' @examples
#' require(lme4)
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' RMSE.merMod(m2)
#' @export
RMSE.merMod <- function(merMod, scale = FALSE){
  stopifnot(inherits(merMod, "lmerMod") ||  inherits(merMod, "blmerMod"))
  # Express RMSE as percentage of dependent variable standard deviation
  dvSD <- sd(merMod@frame[, 1])
  RMSE <- sqrt(mean(residuals(merMod)^2))
  if(scale == TRUE){
    return(RMSE/dvSD)
  } else{
    return(RMSE)
  }
}
