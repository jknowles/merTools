#' Extract theta parameters from a merMod model
#'
#' @param model a valide merMod object
#'
#' @return a vector of the theta parameters from a merMod
#' @export
#'
thetaExtract <- function(model){
  return(model@theta)
}

#' Bootstrap a subset of an lme4 model
#'
#' @param model a valid merMod object
#' @param n the number of rows to sample from the original data
#' in the merMod object
#' @param FUN the function to apply to each bootstrapped model
#' @param B the number of bootstrap replicates
#'
#' @return a matrix of parameters extracted from each of the B replications
#' @details This function allows users to estimate parameters of a
#' large merMod object using bootstraps on a subset of the data.
#' @export
subBoot <- function(model, n, FUN, B){
  resultMat <- matrix(FUN(model), nrow = 1)
  tmp <- matrix(data=NA, nrow=B, ncol=ncol(resultMat))
  resultMat <- rbind(resultMat, tmp); rm(tmp)
    for(i in 1:B){
      rows <- as.numeric(row.names(model@frame))
      mysamp <- as.character(sample(rows,  n, replace=TRUE))
      # http://proteo.me.uk/2013/12/fast-subset-selection-by-row-name-in-r/
      newdata <- model@frame[match(mysamp, rownames(model@frame)),]
      # Only for lmerMod
      tmpMod <- lmer(formula(model), data = newdata)
      resultMat[i + 1, ] <- FUN(tmpMod)
    }
  resultMat <- as.data.frame(resultMat)
  resultMat$sample <- c("original", 1:B)
  return(resultMat)
}


