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

##' Random Effects formula only
##' @keywords internal
reOnly <- function(f,response=FALSE) {
  response <- if (response && length(f)==3) f[[2]] else NULL
  reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
              response=response)
}

safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)

#' Build model matrix
#' @description a function to create a model matrix with all predictor terms in
#' both the group level and fixed effect level
#' @param model a merMod object from lme4
#' @param newdata a data frame to construct the matrix from
#' @param character which matrix to return,default is full matrix with fixed and
#' random terms, other options are "fixed" and "random"
#' @import lme4
#' @keywords internal
buildModelMatrix <- function(model, newdata, which = "full"){
  X <- getME(model, "X")
  X.col.dropped <- attr(X, "col.dropped")
  ## modified from predict.glm ...
  if (is.null(newdata)) {
    newdata <- model@frame
  }
  RHS <- formula(substitute(~R,
                    list(R= RHSForm(formula(model,fixed.only=TRUE)))))
  Terms <- terms(model,fixed.only=TRUE)
  mf <- model.frame(model, fixed.only=TRUE)
  isFac <- vapply(mf, is.factor, FUN.VALUE=TRUE)
  isFac[attr(Terms,"response")] <- FALSE
  orig_levs <- if (length(isFac)==0) NULL else lapply(mf[isFac],levels)
  mfnew <- model.frame(delete.response(Terms),
                         newdata,
                         na.action="na.pass",
                         xlev=orig_levs)
  X <- model.matrix(RHS, data=mfnew,
                      contrasts.arg=attr(X,"contrasts"))
  offset <- 0 # rep(0, nrow(X))
  tt <- terms(model)
  if (!is.null(off.num <- attr(tt, "offset"))) {
    for (i in off.num)
      offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
  }
  fit.na.action <- attr(mfnew,"na.action")
  if(is.numeric(X.col.dropped) && length(X.col.dropped) > 0)
    X <- X[, -X.col.dropped, drop=FALSE]
  re.form <- reOnly(formula(model)) # RE formula only
  newRE <- mkNewReTrms(object = model,
                         newdata = newdata, re.form, na.action="na.pass",
                         allow.new.levels=TRUE)
  reMat <- t(as.matrix(newRE$Zt))
  reMat <- as.matrix(reMat)
  colnames(reMat) <- rownames(newRE$Zt)
  mm <- cbind(X, reMat)
  if(which == "full"){
    return(mm)
  } else if(which == "fixed"){
    return(X)
  } else if(which == "random"){
    return(reMat)
  }

}

#' Calculate the intraclass correlation using mixed effect models
#'
#' @param outcome a character representing the variable of the outcome
#' @param group a character representing the name of the grouping term
#' @param data a data.frame
#' @param subset an optional subset
#'
#' @return a numeric for the intraclass correlation
#' @export
#' @import lme4
#' @examples
#' data(sleepstudy)
#' ICC(outcome = "Reaction", group = "Subject", data = sleepstudy)
ICC <- function(outcome, group, data, subset=NULL){
  fm1 <- as.formula(paste(outcome, "~", "1 + (1|", group, ")"))
  if(length(table(data[, outcome])) == 2){
    nullmod <- glmer(fm1, data = data, subset = subset, family = 'binomial')
  } else {
    nullmod <- lmer(fm1, data = data, subset = subset)
  }
  between <- as.numeric(attr(VarCorr(nullmod)[[1]], "stddev"))
  within <- sigma(nullmod)
  ICC <- between / (within + between)
  return(ICC)
}


#' Utility function to make RE terms objects
#' @import lme4
#' @keywords internal
mkNewReTrms <- function(object, newdata, re.form=NULL, na.action=na.pass,
                        allow.new.levels=FALSE)
{
  if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    mfnew <- model.frame(delete.response(terms(object, fixed.only=TRUE)),
                         newdata, na.action=na.action)
    if(packageVersion("lme4") < "1.1.9"){
      old <- TRUE
    } else{
      old <- FALSE
    }

    if (old) {
      rfd <- na.action(newdata)
      if (is.null(attr(rfd,"na.action")))
        attr(rfd,"na.action") <- na.action
    } else {
      newdata.NA <- newdata
      if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
        newdata.NA <- newdata.NA[-fixed.na.action,]
      }
      tt <- delete.response(terms(object,random.only=TRUE))
      ## need to let NAs in RE components go through -- they're handled downstream
      rfd <- model.frame(tt,newdata.NA,na.action=na.pass)
      if (!is.null(fixed.na.action))
        attr(rfd,"na.action") <- fixed.na.action
    }
  }
  if (inherits(re.form, "formula")) {
    ## DROP values with NAs in fixed effects
    if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
      newdata <- newdata[-fit.na.action,]
    }
    ## note: mkReTrms automatically *drops* unused levels
    ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
    ## update Lambdat (ugh, better way to do this?)
    ReTrms <- within(ReTrms,Lambdat@x <- unname(getME(object,"theta")[Lind]))
    if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
      stop("NAs are not allowed in prediction data",
           " for grouping variables unless allow.new.levels is TRUE")
    ns.re <- names(re <- ranef(object))
    nRnms <- names(Rcnms <- ReTrms$cnms)
    if (!all(nRnms %in% ns.re))
      stop("grouping factors specified in re.form that were not present in original model")
    new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
    ## fill in/delete levels as appropriate
    re_x <- Map(function(r,n) levelfun(r,n,allow.new.levels=allow.new.levels),
                re[names(new_levels)], new_levels)
    re_new <- lapply(seq_along(nRnms), function(i) {
      rname <- nRnms[i]
      if (!all(Rcnms[[i]] %in% names(re[[rname]])))
        stop("random effects specified in re.form that were not present in original model")
      re_x[[rname]][,Rcnms[[i]]]
    })
    re_new <- unlist(lapply(re_new, t))
  }
  Zt <- ReTrms$Zt
  attr(Zt, "na.action") <- attr(re_new, "na.action") <- attr(mfnew, "na.action")
  list(Zt=Zt, b=re_new, Lambdat = ReTrms$Lambdat)
}

#' Parse merMod formulas
#' @keywords internal
RHSForm <- function(form,as.form=FALSE) {
  rhsf <- form[[length(form)]]
  if (as.form) reformulate(deparse(rhsf)) else rhsf
}


#' Parse merMod levels
#' @keywords internal
levelfun <- function(x,nl.n,allow.new.levels=FALSE) {
  if (!all(nl.n %in% rownames(x))) {
    if (!allow.new.levels) stop("new levels detected in newdata")
      newx <- as.data.frame(matrix(0, nrow=length(nl.n), ncol=ncol(x),
                                 dimnames=list(nl.n, names(x))))
      newx[rownames(x),] <- x
      x <- newx
  }
  if (!all(r.inn <- rownames(x) %in% nl.n)) {
    x <- x[r.inn,,drop=FALSE]
  }
  return(x)
}
