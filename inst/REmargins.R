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
REmargins <- function(merMod, newdata, groupFctr = NULL, term = NULL, breaks = 3,
                      .parallel = TRUE, ...){
  if(missing(groupFctr)){
    groupFctr <- names(ranef(merMod))[1]
  }

  # This is a rough way to break the ER distribution into quantiles
  brks <- ceiling(seq(1, 100, by = 100/breaks))
  if(! 99 %in% brks) {
    brks <- c(brks, 99)
  }



  # Generate the expected rank distribution
  ER_DF <- expectedRank(merMod, groupFctr = groupFctr)
  ER_DF <- ER_DF[!duplicated(ER_DF[, c("groupFctr", "term", "pctER")]), ]
  # Keep only factor levels that have effects at the margins
  # Need to match closest value here
  # Find N closest values
  # Drop duplicates
  zz <- abs(diff(floor(c(brks, ER_DF$pctER))))
  ndx <- order(zz, decreasing = TRUE)[1:5]

  ER_DF <- ER_DF[ndx, ]
  # Question: Should we instead prefer the most precise effect at each quantile?
  # Find all the levels as a list
  lvls <- as.data.frame(unique(merMod@frame[, groupFctr]))
  names(lvls) <- groupFctr
  # Set up a list to step through
  keep <- vector(mode = "list", length(names(lvls)))
  names(keep) <- names(lvls)
  # For each factor level, step through and keep all combinations of terms where
  # that factor level matches one remaining in ER_DF
  for(i in names(lvls)) {
    keep[i] <- lvls[lvls[, i] %in% ER_DF$groupLevel[ER_DF$groupFctr == i], , drop=FALSE]
  }
  # Keep everything that remains
  lvls <- dplyr::bind_rows(keep)
  # Clean up
  rm(keep, ER_DF, brks)

  # TODO - order case by magnitude of ER


  # Get ready to expand the data
  zed <- as.data.frame(lapply(newdata, rep, nrow(lvls)))
  zed[, groupFctr] <- lvls[rep(seq_len(nrow(lvls)), nrow(newdata)), ]
  zed[, "case"] <- rep(seq(1, nrow(newdata)), times = nrow(lvls))

  # Maybe strongly recommend parallel here?
  if ( .parallel & requireNamespace("foreach", quietly=TRUE)) {
    # TODO use future here

      merTools:::setup_parallel()

      fe_call <- as.call(c(list(quote(foreach::foreach), i = unique(zed$case),
                                .inorder = FALSE)))
      fe <- eval(fe_call)
      keep <- foreach::`%dopar%`(fe, {
        out <- predictInterval(merMod, newdata = zed[zed$case == i, ],
                               which = "all")
        out_w <- stats::reshape(out, direction = "wide",
                                idvar = "obs", timevar = "effect",
                                v.names = c("fit", "upr", "lwr"), sep = "_")
        out_w$case <- i
        out_w
      }
      )

      out_w <- dplyr::bind_rows(keep); rm(keep)
    } else if ( .parallel & !requireNamespace("foreach", quietly=TRUE)) {
      warning("foreach package is unavailable, parallel computing not available")
    } else {
    keep <- vector(mode = "list", nrow(newdata))
    for (i in unique(zed$case)) {
      out <- predictInterval(merMod, newdata = zed[zed$case == i, ],
                             which = "all")
      out_w <- stats::reshape(out, direction = "wide",
                       idvar = "obs", timevar = "effect",
                       v.names = c("fit", "upr", "lwr"), sep = "_")
      out_w$case <- i
      keep[[i]] <- out_w
    }
    out_w <- dplyr::bind_rows(keep); rm(keep)
  }

  # Case is the number of rows in newdata
  # obs is the variance among the selected random effects to marginalize over
  # So we want to collapse by case if we can

  return(out_w)
}
