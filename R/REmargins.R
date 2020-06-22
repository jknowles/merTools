#' Calculate the predicted value for each observation across the distribution
#' of the random effect terms.
#'
#' \code{REmargins} calculates the average predicted value for each row of a
#' new data frame across the distribution of \code{\link{expectedRank}} for a
#' merMod object. This allows the user to make meaningful comparisons about the
#' influence of random effect terms on the scale of the response variable,
#' for user-defined inputs, and accounting for the variability in grouping terms.
#'
#' The function simulates the
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
#' in the newdata element:
#'   \describe{
#'     \item{...}{The columns of the original data taken from \code{newdata}}
#'     \item{case}{The row number of the observation from newdata. Each row in newdata will be
#'     repeated for all unique levels of the grouping_var, term, and breaks.}
#'     \item{grouping_var}{The grouping variable the random effect is being marginalized over.}
#'     \item{term}{The term for the grouping variable the random effect is being marginalized over.}
#'     \item{breaks}{The ntile of the effect size for \code{grouping_var} and \code{term}}
#'     \item{original_group_level}{The original grouping value for this \code{case}}
#'     \item{fit_combined}{The predicted value from \code{predictInterval} for this case simulated
#'     at the Nth ntile of the expected rank distribution of \code{grouping_var} and \code{term}}
#'     \item{upr_combined}{The upper bound of the predicted value.}
#'     \item{lwr_combined}{The lower bound of the predicted value.}
#'     \item{fit_XX}{For each grouping term in newdata the predicted value is decomposed into its
#'     fit components via predictInterval and these are all returned here}
#'     \item{upr_XX}{The upper bound for the effect of each grouping term}
#'     \item{lwr_XX}{The lower bound for the effect of each grouping term}
#'     \item{fit_fixed}{The predicted fit with all the grouping terms set to 0 (average)}
#'     \item{upr_fixed}{The upper bound fit with all the grouping terms set to 0 (average)}
#'     \item{lwr_fixed}{The lower bound fit with all the grouping terms set to 0 (average)}
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
#' \donttest{
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:10,])
#'
#' # You can also pass additional arguments to predictInterval through REimpact
#'  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
#'  margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("s"),
#'                         breaks = 4)
#'  margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("d"),
#'                          breaks = 3)
#' }
#' @export
REmargins <- function(merMod, newdata = NULL, groupFctr = NULL, term = NULL, breaks = 4,
                      .parallel = FALSE, ...){
  # Validate inputs
  if (is.null(groupFctr)) {
    # If the user doesn't tell us which term to use, we take the first term
    groupFctr <- names(ranef(merMod))[1]
  }
  if (is.null(newdata)) {
    # If the user doesn't give us data, we take the whole dataset
    # TODO - how does performance scale to a large number of observations?
    newdata <- merMod@frame
  }

  # If the user doesn't tell us what term to use, we take all the terms
  if (is.null(term)) {
    term <- names(ranef(merMod)[[groupFctr]])
    # Sub out intercept
    term[term == "(Intercept)"] <- "Intercept"
  }

  # This is a rough way to break the ER distribution into quantiles
  brks <- ceiling(seq(1, 100, by = 100/breaks))
  # Fallback so we always take a 99th percentile effect (for the maximum)
  if (!99 %in% brks) {
    brks <- c(brks, 99)
  }
  # Inputs are validated - now we get the effect distribution
  # Generate the expected rank distribution
  ER_DF <- expectedRank(merMod, groupFctr = groupFctr, term = term)
  # With many effects there is a lot of duplication - drop duplicated pctER
  ER_DF <- ER_DF[!duplicated(ER_DF[, c("groupFctr", "term", "pctER")]), ]

  # Now we create a data frame to capture the factor levels of each groupFctr that
  # correspond to the right break in the expectedRank distribution of the random
  # effect grouping factor and term
  par_df <- expand.grid("grouping_var" = groupFctr, "term" = term, "breaks" = 1:breaks)

  # Keep only factor levels that have effects at the margins
  # Need to match closest value here
  # Find N closest values
  # Drop duplicates
    # For each combination build an index of candidate rows/effect levels
  # Then choose the level that has the most precise estimate within a
  # tolerance of the effect size
  for (trm in term) {
    for (i in seq_along(brks)) {
      # Compute each terms distance from the break
      rank_dist <- abs(ER_DF$pctER[ER_DF$term == trm] - brks[i])
      # Get the index for the rank that minimizes the distance
      # TODO - how to break ties here?
      tmp <- which(rank_dist %in% rank_dist[order(rank_dist)][1])
      # Store the result in the par_df object
      par_df$groupLevel[par_df$breaks == i & par_df$term == trm] <- ER_DF$groupLevel[tmp]
    }
  }
  # Get ready to expand the data
  sim_data <- as.data.frame(lapply(newdata, rep, each = nrow(par_df)))
  # sim_data now repeats each row of newdata by the number of rows in par_df
  # case labels the rows with an integer for later mapping
  sim_data$case <- rep(1:nrow(newdata), each = nrow(par_df))
  sim_data <- cbind(sim_data, par_df)
  sim_data$original_group_level <- sim_data[, groupFctr]
  sim_data[, groupFctr] <- sim_data$groupLevel
  sim_data$groupLevel <- NULL
  #

  # Maybe strongly recommend parallel here?
  if (.parallel & requireNamespace("foreach", quietly = TRUE)) {
    # TODO use future here
      setup_parallel()
        out <- predictInterval(merMod, newdata = sim_data,
                               which = "all", ...)
        out_w <- stats::reshape(out, direction = "wide",
                                idvar = "obs", timevar = "effect",
                                v.names = c("fit", "upr", "lwr"), sep = "_")
        out_w$obs <- NULL
        sim_data <- cbind(sim_data, out_w)
    } else if ( .parallel & !requireNamespace("foreach", quietly = TRUE)) {
      warning("foreach package is unavailable, parallel computing not available")
    } else {
      out <- predictInterval(merMod, newdata = sim_data,
                             which = "all", ...)
      out_w <- stats::reshape(out, direction = "wide",
                       idvar = "obs", timevar = "effect",
                       v.names = c("fit", "upr", "lwr"), sep = "_")
      out_w$obs <- NULL
      sim_data <- cbind(sim_data, out_w)
    }
  # Case is the number of the row in newdata
  # obs is the variance among the selected random effects to marginalize over
  # So we want to collapse by case if we can

  return(sim_data)
}
