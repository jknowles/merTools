#' Calculate the weighted mean of fitted values for various combinations of
#' random effect terms.
#'
#' \code{REmargins} calculates the average predicted value for each row of a
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
#'   by \code{ranef}. If the length is > 1 then the combined effect of all
#'   listed groups will calculated and marginalized over co-occurences of those
#'   groups if desired.
#'
#' @param term The name of the random coefficient of interest. This is the
#'   variable to the left of the pipe, \code{|}, in the [g]lmer formula. Partial
#'   matching is attempted on the intercept term so the following character
#'   strings will all return rankings based on the intercept (\emph{provided that
#'   they do not match the name of another random coefficient for that factor}):
#'   \code{c("(Intercept)", "Int", "intercep", ...)}.
#'
#' @param breaks an integer representing the number of bins to divide the group
#' effects into, the default is 3.
#' @param .parallel, logical should parallel computation be used, default is TRUE
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
#' @importFrom stats reshape
#' @examples
#' #For a one-level random intercept model
#' \donttest{
#' # You can also pass additional arguments to predictInterval through REimpact
#'  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
#'  margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("s"),
#'                         breaks = 4)
#'  margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("d"),
#'                          breaks = 3)
#' }
#' @export
REmargins <- function(merMod, newdata = NULL, groupFctr = NULL, term = NULL, breaks = 3,
                      .parallel = FALSE, ...){
  if(is.null(groupFctr)){
    groupFctr <- names(ranef(merMod))[1]
  }
  if (is.null(newdata)) {
    newdata <- merMod@frame
  }

  # This is a rough way to break the ER distribution into quantiles
  brks <- ceiling(seq(1, 100, by = 100/breaks))
  if(! 99 %in% brks) {
    brks <- c(brks, 99)
  }
  # Generate the expected rank distribution
  ER_DF <- expectedRank(merMod, groupFctr = groupFctr)
  ER_DF <- ER_DF[!duplicated(ER_DF[, c("groupFctr", "term", "pctER")]), ]


  if (is.null(term)) {
    term <- unique(ER_DF$term)
  }

  out_df <- expand.grid("grouping_var" = groupFctr, "term" = term, "breaks" = 1:breaks)
  out_df$groupLevel <- NA


  # Keep only factor levels that have effects at the margins
  # Need to match closest value here
  # Find N closest values
  # Drop duplicates
    # For each combination build an index of candidate rows/effect levels
  # Then choose the level that has the most precise estimate within a
  # tolerance of the effect size


  for (trm in term) {
    for (i in seq_along(brks)) {

      zz <- abs(ER_DF$pctER[ER_DF$term == trm] - brks[i])
      tmp <- which(zz %in% zz[order(zz)][1])
      out_df$groupLevel[out_df$breaks == i & out_df$term == trm] <- ER_DF$groupLevel[tmp]
    }
  }
  # Need to repeat
  # TODO - Figure out how to keep the effects seperate from each term
  # TODO - order case by magnitude of ER

  # Get ready to expand the data
  zed <- as.data.frame(lapply(newdata, rep, each = nrow(out_df)))
  zed$case <- rep(1:nrow(newdata), each = nrow(out_df))
  zed <- cbind(zed, out_df)
  zed$original_group_level <- zed[, groupFctr]
  zed[, groupFctr] <- zed$groupLevel
  zed$groupLevel <- NULL
  #
  # zed[, "case"] <- rep(seq(1, nrow(newdata)), times = nrow(lvls))

  # Maybe strongly recommend parallel here?
  if ( .parallel & requireNamespace("foreach", quietly=TRUE)) {
    # TODO use future here
      setup_parallel()
        out <- predictInterval(merMod, newdata = zed,
                               which = "all")
        out_w <- stats::reshape(out, direction = "wide",
                                idvar = "obs", timevar = "effect",
                                v.names = c("fit", "upr", "lwr"), sep = "_")
        out_w$obs <- NULL
        zed <- cbind(zed, out_w)
    } else if ( .parallel & !requireNamespace("foreach", quietly=TRUE)) {
      warning("foreach package is unavailable, parallel computing not available")
    } else {
      out <- predictInterval(merMod, newdata = zed,
                             which = "all")
      out_w <- stats::reshape(out, direction = "wide",
                       idvar = "obs", timevar = "effect",
                       v.names = c("fit", "upr", "lwr"), sep = "_")
      out_w$obs <- NULL
      zed <- cbind(zed, out_w)
    }
  # Case is the number of the row in newdata
  # obs is the variance among the selected random effects to marginalize over
  # So we want to collapse by case if we can

  return(zed)
}
