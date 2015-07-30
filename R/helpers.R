#Helpers

# Function to take only rows that form distinct levels of factors
# Need to figure out how to build a model matrix better.
trimModelFrame <- function(data){
  # Identify numerics
  nums <- sapply(data, is.numeric)
  vars <- names(nums[!nums == TRUE])
  dataList <- vector(mode = "list", length = length(vars))
  names(dataList) <- vars
    for(i in vars){
      dataList[[i]] <- data[!duplicated(data[, i]), ,drop=FALSE]
    }
    newdat <- do.call(rbind, dataList)
    newdat <- newdat[!duplicated(newdat),]
    return(newdat)
}

# FROM LME4
residDF.merMod <- function(object) {
  npar <- length(object@beta) + length(object@theta) +
    object@devcomp[["dims"]][["useSc"]]
  nobs <- nrow(object@frame)
  ## TODO: how do we feel about counting the scale parameter ???
  return(nobs - npar)
}


# from ARM as.matrix.VarCorr
easyVarCorr <- function (varc, useScale, digits){
  # VarCorr function for lmer objects, altered as follows:
  #   1.  specify rounding
  #   2.  print statement at end is removed
  #   3.  reMat is returned
  #   4.  last line kept in reMat even when there's no error term
  sc <- attr(varc, "sc")[[1]]
  if(is.na(sc)) sc <- 1
  #                  recorr <- lapply(varc, function(el) el@factors$correlation)
  recorr <- lapply(varc, function(el) attr(el, "correlation"))
  #reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
  reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
  reLens <- unlist(c(lapply(reStdDev, length)))
  reMat <- array('', c(sum(reLens), 4),
                 list(rep('', sum(reLens)),
                      c("Groups", "Name", "Variance", "Std.Dev.")))
  reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
  reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
  #                  reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
  #                  reMat[,4] <- format(unlist(reStdDev), digits = digits)
  reMat[,3] <- fround(unlist(reStdDev)^2, digits)
  reMat[,4] <- fround(unlist(reStdDev), digits)
  if (any(reLens > 1)) {
    maxlen <- max(reLens)
    corr <-
      do.call("rbind",
              lapply(recorr,
                     function(x, maxlen) {
                       x <- as(x, "matrix")
                       #   cc <- format(round(x, 3), nsmall = 3)
                       cc <- fround (x, digits)
                       cc[!lower.tri(cc)] <- ""
                       nr <- dim(cc)[1]
                       if (nr >= maxlen) return(cc)
                       cbind(cc, matrix("", nr, maxlen-nr))
                     }, maxlen))
    colnames(corr) <- c("Corr", rep("", maxlen - 1))
    reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
  }
  #    if (!useScale) reMat <- reMat[-nrow(reMat),]
  if (useScale<0) reMat[nrow(reMat),] <- c ("No residual sd", rep("",ncol(reMat)-1))
  return (reMat)
}

reTermCount <- function(model){
  sum(unlist(lapply(as.list(VarCorr(model)), function(x) sqrt(length(x)))))
}

reTermNames <- function(model){
  tmp <- NA
  for(i in 1:length(names(ngrps(model)))){
    cons <- names(ngrps(model))[i]
    vars <- paste(cons, unlist(dimnames(VarCorr(model)[[i]])[1]), sep = "-")
    tmp <- c(tmp, vars)
  }
  tmp <- na.omit(tmp)
  tmp <- t(as.data.frame(strsplit(tmp, "-")))
  row.names(tmp) <- NULL
  colnames(tmp) <- c("group", "effect")
  tmp <- as.data.frame(tmp)
  tmp$group <- as.character(tmp$group)
  tmp$effect <- as.character(tmp$effect)
  return(tmp)
}

#' Clean formula
#' @description a function to modify the formula for a merMod object to create
#' a model matrix with all predictor terms in both the group level and fixed
#' effect level
#' @param model a merMod object from lme4
#' @keywords internal
formulaBuild <- function(model){
  slopeFX <- setdiff(all.vars(model@call$formula), names(ngrps(model)))
  missVar <- setdiff(slopeFX, all.vars(nobars(model@call$formula)))
  newForm <- nobars(model@call$formula)
  if(length(missVar > 0)){
    newForm <- paste(Reduce(paste, deparse(newForm)), paste(missVar, collapse = " +"), sep = " + ")
  }
  newForm <- as.formula(newForm)
  return(newForm)
}




