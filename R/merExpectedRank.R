#' Calculate the expected rank of random coefficients that account for
#' uncertainty.
#'
#' \code{expectedRank} calculates the expected rank and the percentile expected
#' rank of any random term in a merMod object.  A simple ranking of the estimated
#' random effects (as produced by \code{\link[lme4]{ranef}}) is not satisfactory
#' because it ignores any amount of uncertainty.
#'
#' Inspired by Lingsma et al. (2010, see also Laird and Louis 1989),
#' expectedRank sums the probability that each level of the grouping factor is
#' greater than every other level of the grouping factor, similar to a
#' two-sample t-test.
#'
#' The formula for the expected rank is:
#' \deqn{ExpectedRank_i = 1 + \sum \phi((\theta_i - \theta_k) / \sqrt(var(\theta_i)+var(\theta_k))}
#' where \eqn{\phi} is the standard normal distribution function, \eqn{\theta}
#' is the estimated random effect and \eqn{var(\theta)} is the posterior
#' variance of the estimated random effect. We add one to the sum so that the
#' minimum rank is one instead of zero so that in the case where there is no
#' overlap between the variances of the random effects (or if the variances are
#' zero), the expected rank equals the actual rank.  The ranks are ordered such
#' that the winners have ranks that are greater than the losers.
#'
#' The formula for the percentile expected rank is:
#' \deqn{100 * (ExpectedRank_i - 0.5) / N_grps}
#' where \eqn{N_grps} is the number of grouping factor levels. The percentile
#' expected rank can be interpreted as the fraction of levels that score at or
#' below the given level.
#'
#' NOTE: \code{expectedRank} will only work under conditions that \code{lme4::ranef}
#' will work. One current example of when this is \emph{not} the case is for
#' models when there are multiple terms specified per factor (e.g. uncorrelated random
#' coefficients for the same term, e.g.
#' \code{lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data = sleepstudy)})
#'
#' @param merMod An object of class merMod
#'
#' @param factor The name of the grouping factor over which the random
#'   coefficient of interest varies.  This is the variable to the right of the
#'   pipe, \code{|}, in the [g]lmer formula. This parameter is optional if only
#'   a single grouping factor is included in the model, but required if there
#'   are two or more.
#'
#' @param term The name of the random coefficient of interest. This is the
#'   variable to the left of the pipe, \code{|}, in the [g]lmer formula. Partial
#'   matching is attempted on the intercept term so the following character
#'   strings will all return rankings based on the intercept (\emph{provided that
#'   they do not match the name of another random coefficient for that factor}):
#'   \code{c("(Intercept)", "Int", "intercep", ...)}.
#'
#' @return A data.frame with the original grouping factor (converted to a
#'   character) and the following four columns:
#'   \describe{
#'     \item{theta}{The estimated random effect (from
#'                  \code{lme4::ranef(, condVar=TRUE)}).}
#'     \item{varTheta}{The posterior variance of the estimate random effect
#'                      (from \code{lme4::ranef(, condVar=TRUE)}).}
#'     \item{ER}{The expected rank.}
#'     \item{pctER}{The percentile expected rank.}
#'   }
#'
#' @references
#' Laird NM and Louis TA. Empirical Bayes Ranking Methods. \emph{Journal of
#' Education Statistics}. 1989;14(1)29-46. Available at
#' \url{http://www.jstor.org/stable/1164724}.
#'
#'
#' Lingsma HF, Steyerberg EW, Eijkemans MJC, et al. Comparing and
#' ranking hospitals based on outcome: results from The Netherlands Stroke Survey.
#' \emph{QJM: An International Journal of Medicine}. 2010;103(2):99-108.
#' doi:10.1093/qjmed/hcp169
#'
#' @examples
#' #For a one-level random intercept model
#' require(lme4)
#' m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
#' (m1.er <- expectedRank(m1))
#'
#' #For a one-level random intercept model with multiple random terms
#' require(lme4)
#' m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' #ranked by the random slope on Days
#' (m2.er1 <- expectedRank(m2, term="Days"))
#' #ranked by the random intercept
#' (m2.er2 <- expectedRank(m2, term="int"))
#'
#' \dontrun{
#' #For a two-level model with random intercepts
#' require(lme4)
#' m3 <- lmer(y ~ service * dept + (1|s) + (1|d), InstEval)
#' #Ranked by the random intercept on 's'
#' (m3.er1 <- expectedRank(m3, factor="s", term="Intercept"))
#' }
#' @export
expectedRank <- function(merMod, factor=NULL, term=NULL) {
  #Count random terms in merMod
  n.rfx <- lme4::getME(merMod, "k")
  n.rfac <- lme4::getME(merMod, "n_rfac")

  rfx <- lme4::ranef(merMod, condVar=TRUE)

  #Take care of factors
  if (n.rfac == 1 & is.null(factor)) {
    factor <- names(rfx)
  }
  else if (n.rfac > 1 & is.null(factor)) {
    stop("Must specify which grouping factor when there are more than one")
  }

  factor.idx <- grep(factor, names(rfx), fixed=TRUE)

  #Grab row names, number of columns
  rfx.names <- rownames(rfx[[factor.idx]])
  n.grps <- length(rfx.names)
  n.terms <- length(rfx[[factor.idx]])

  #Take care of term
  if (n.terms == 1 & is.null(term)) {
    term <- names(rfx[[factor.idx]])
  }
  else if (n.terms > 1 & is.null(term)) {
    stop("Must specify which random coefficient when there are more than one per selected grouping factor")
  }

  if (grepl("[iI]nt[a-z]*", term) && is.na(match(term, names(rfx[[factor.idx]])))) {
    term <- "(Intercept)"
  }

  term.idx <- grep(term, names(rfx[[factor.idx]]), fixed=TRUE)

  #Grab theta and var.theta
  theta <- rfx[[factor.idx]][,term.idx]
  var.theta <- attr(rfx[[factor.idx]], which="postVar")[term.idx, term.idx, 1:n.grps]

  #Calculate Expected Rank which is the sum of the probabilities that group i is greater than all
  #other groups j (assuming normal distribution of random effects)
  ER <- pctER <- rep(NA, n.grps)
  for (i in 1:n.grps) {
    ER[i] <- 1 + sum(pnorm((theta[i]-theta[-i]) / sqrt(var.theta[i] + var.theta[-i])))
  }

  #Calculated percentile expected rank ... the version of the formula I am using is
  #the percentage of groups that are ranked **equal to or less than** the selected
  #group ... if we just wanted percentage ranked less than then remove the 0.5
  pctER <- round(100 * (ER - 0.5) / n.grps)

  #Close out and return in order of best to worst
  out <- data.frame(rfx.names, theta, var.theta, ER, pctER, stringsAsFactors=FALSE)
  names(out) <- c(factor, "theta", "varTheta", "ER", "pctER")
  out <- out[order(out$ER, decreasing=TRUE),]
  row.names(out) <- 1:n.grps
  return(out)
}
