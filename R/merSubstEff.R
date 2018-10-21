#' Calculate the weighted mean of fitted values for various levels of
#' random effect terms.
#'
#' \code{REimpact} calculates the average predicted value for each row of a
#' new data frame across the distribution of \code{\link{expectedRank}} for a
#' merMod object. This allows the user to make meaningful comparisons about the
#' influence of random effect terms on the scale of the response variable,
#' for user-defined inputs, and accounting for the variability in grouping terms.
#'
#' The function predicts the response at every level in the random effect term
#' specified by the user. Then, the expected rank of each group level is binned
#' to the number of bins specified by the user. Finally, a weighted mean of the
#' fitted value for all observations in each bin of the expected ranks is
#' calculated using the inverse of the variance as the weight -- so that less
#' precise estimates are downweighted in the calculation of the mean for the bin.
#' Finally, a standard error for the bin mean is calculated.
#'
#' @param merMod An object of class merMod
#'
#' @param newdata a data frame of observations to calculate group-level differences
#' for
#'
#' @param groupFctr The name of the grouping factor over which the random
#'   coefficient of interest varies.  This is the variable to the right of the
#'   pipe, \code{|}, in the [g]lmer formula. This parameter is optional, if not
#'   specified, it will perform the calculation for the first effect listed
#'   by \code{ranef}.
#'
#' @param term The name of the random coefficient of interest. This is the
#'   variable to the left of the pipe, \code{|}, in the [g]lmer formula. Partial
#'   matching is attempted on the intercept term so the following character
#'   strings will all return rankings based on the intercept (\emph{provided that
#'   they do not match the name of another random coefficient for that factor}):
#'   \code{c("(Intercept)", "Int", "intercep", ...)}.
#'
#' @param breaks an integer representing the number of bins to divide the group
#' effects into, the default is 3; alternatively it can specify breaks from 0-100
#' for how to cut the expected rank distribution
#'
#' @param ... additional arguments to pass to \code{\link{predictInterval}}
#'
#' @return A data.frame with all unique combinations of the number of cases, rows
#' in the newdata element, and number of bins:
#'   \describe{
#'     \item{case}{The row number of the observation from newdata.}
#'     \item{bin}{The ranking bin for the expected rank, the higher the bin number,
#'     the greater the expected rank of the groups in that bin.}
#'     \item{AvgFitWght}{The weighted mean of the fitted values for case i in bin k}
#'     \item{AvgFitWghtSE}{The standard deviation of the mean of the fitted values
#'     for case i in bin k.}
#'     \item{nobs}{The number of group effects contained in that bin.}
#'   }
#'
#' @details This function uses the formula for variance of a weighted mean
#' recommended by Cochran (1977).
#'
#' @references
#' Gatz, DF and Smith, L. The Standard Error of a Weighted Mean Concentration.
#' I. Bootstrapping vs other methods. \emph{Atmospheric Environment}.
#' 1995;11(2)1185-1193. Available at
#' \url{http://www.sciencedirect.com/science/article/pii/135223109400210C}
#'
#' Cochran, WG. 1977. Sampling Techniques (3rd Edition). Wiley, New York.
#'
#' @seealso \code{\link{expectedRank}}, \code{\link{predictInterval}}
#'
#' @examples
#' #For a one-level random intercept model
#' m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
#' m1.er <- REimpact(m1, newdata = sleepstudy[1, ], breaks = 2)
#' #For a one-level random intercept model with multiple random terms
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' #ranked by the random slope on Days
#' m2.er1 <- REimpact(m2,  newdata = sleepstudy[1, ],
#'            groupFctr = "Subject", term="Days")
#' #ranked by the random intercept
#' m2.er2 <- REimpact(m2, newdata = sleepstudy[1, ],
#'              groupFctr = "Subject", term="int")
#' \donttest{
#' # You can also pass additional arguments to predictInterval through REimpact
#' g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
#' zed <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", n.sims = 50,
#'                 include.resid.var = TRUE)
#' zed2 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "s", n.sims = 50,
#'                  include.resid.var = TRUE)
#' zed3 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", breaks = 5,
#                 n.sims = 50, include.resid.var = TRUE)
#' }
#' @export
REimpact <- function(merMod, newdata, groupFctr=NULL, term = NULL, breaks = 3, ...){
  if(missing(groupFctr)){
    groupFctr <- names(ranef(merMod))[1]
  }
  lvls <- unique(merMod@frame[, groupFctr])
  zed <- as.data.frame(lapply(newdata, rep, length(lvls)))
  zed[, groupFctr] <- rep(lvls, each = nrow(newdata))
  zed[, "case"] <- rep(seq(1, nrow(newdata)), times = length(lvls))
  outs1 <- cbind(zed, predictInterval(merMod, newdata = zed, ...))
  outs1$var <- outs1$upr - outs1$lwr
  outs1$lwr <- NULL; outs1$upr <- NULL
  ranks <- expectedRank(merMod, groupFctr = groupFctr, term = term)
  ranks <- ranks[, c(2, 7)]
  outs1 <- merge(ranks, outs1, by.x = "groupLevel", by.y = groupFctr); rm(ranks)

  weighted.var.se <- function(x, w, na.rm=FALSE)
    #  Computes the variance of a weighted mean following Cochran 1977 definition
  {
    if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
    n = length(w)
    xWbar = weighted.mean(x,w,na.rm=na.rm)
    wbar = mean(w)
    out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
    return(out)
  }
  # bin pctER somehow
  outs1$bin <- cut(outs1$pctER, breaks = breaks, labels = FALSE,
                   include.lowest = TRUE)
  bySum <- function(x){
    AvgFit <- weighted.mean(x$fit, 1/x$var)
    AvgFitSE <- weighted.var.se(x$fit, 1/x$var)
    nobs <- length(x$fit)
    return(c(AvgFit, AvgFitSE, nobs))
  }
  outs1 <- outs1[order(outs1$case, outs1$bin),]

  wMeans <- by(outs1, INDICES = list(outs1$case, outs1$bin), bySum)
  ids <- expand.grid(unique(outs1$case), unique(outs1$bin))
  wMeans <- cbind(ids, do.call(rbind, wMeans))
  names(wMeans) <- c("case", "bin", "AvgFit", "AvgFitSE", "nobs")
  return(wMeans)
}
