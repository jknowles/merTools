#' Extract theta parameters from a merMod model
#' @description A convenience function that returns the theta parameters for a
#' \code{\link{merMod}} obejct.
#' @param merMod a valide merMod object
#'
#' @return a vector of the covrariance, theta, parameters from a \code{\link{merMod}}
#' @seealso merMod
#' @export
#' @examples
#' (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#' thetaExtract(fm1) #(a numeric vector of the covariance parameters)
thetaExtract <- function(merMod){
  stopifnot(class(merMod) %in% c("lmerMod", "glmerMod", "blmerMod",
                                 "bglmerMod"))
  return(merMod@theta)
}

#' Bootstrap a subset of an lme4 model
#'
#' @param merMod a valid merMod object
#' @param n the number of rows to sample from the original data
#' in the merMod object, by default will resample the entire model frame
#' @param FUN the function to apply to each bootstrapped model
#' @param R the number of bootstrap replicates, default is 100
#' @param seed numeric, optional argument to set seed for simulations
#' @return a data.frame of parameters extracted from each of the R replications.
#' The original values are appended to the top of the matrix.
#' @details This function allows users to estimate parameters of a
#' large merMod object using bootstraps on a subset of the data.
#' @examples
#' \dontrun{
#' (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#' resultMatrix <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20)
#' }
#' @export
subBoot <- function(merMod, n = NULL, FUN, R = 100, seed=NULL){
  if(missing(n))
    n <- nrow(merMod@frame)
  resultMat <- matrix(FUN(merMod), nrow = 1)
  tmp <- matrix(data=NA, nrow=R, ncol=ncol(resultMat))
  resultMat <- rbind(resultMat, tmp); rm(tmp)

  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

    for(i in 1:R){
      rows <- as.numeric(row.names(merMod@frame))
      mysamp <- as.character(sample(rows,  n, replace=TRUE))
      # http://proteo.me.uk/2013/12/fast-subset-selection-by-row-name-in-r/
      newdata <- merMod@frame[match(mysamp, rownames(merMod@frame)),]
      # Only for lmerMod
      tmpMod <- lmer(formula(merMod), data = newdata)
      resultMat[i + 1, ] <- FUN(tmpMod)
    }
  resultMat <- data.frame(param=resultMat)
  resultMat$replicate <- c("original", 1:R)
  return(resultMat)
}


