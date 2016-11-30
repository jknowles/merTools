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
#' @param groupFctr An optional character vector specifying the name(s) the grouping factor(s)
#'   over which the random coefficient of interest varies.  This is the
#'   variable to the right of the pipe, \code{|}, in the [g]lmer formula.
#'   This parameter is optional. If none is specified all terms will be returned.
#'
#' @param term An optional character vector specifying the name(s) of the random coefficient of interest. This is the
#'   variable to the left of the pipe, \code{|}, in the [g]lmer formula. Partial
#'   matching is attempted on the intercept term so the following character
#'   strings will all return rankings based on the intercept (\emph{provided that
#'   they do not match the name of another random coefficient for that factor}):
#'   \code{c("(Intercept)", "Int", "intercep", ...)}.
#'
#' @return A data.frame with the following five columns:
#'   \describe{
#'     \item{groupFctr}{a character representing name of the grouping factor}
#'     \item{groupLevel}{a character representing the level of the grouping factor}
#'     \item{term}{a character representing the formula term for the group}
#'     \item{estimate}{effect estimate from \code{lme4::ranef(, condVar=TRUE)}).}
#'     \item{std.error}{the posterior variance of the estimate random effect
#'                      (from \code{lme4::ranef(, condVar=TRUE)}); named "\code{term}"_var.}
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
#' (m3.er1 <- expectedRank(m3, groupFctr="s", term="Intercept"))
#' }
#' @export
expectedRank <- function(merMod, groupFctr=NULL, term=NULL) {
  #Count random terms in merMod
  n.rfx <- lme4::getME(merMod, "k")
  n.rfac <- lme4::getME(merMod, "n_rfac")
  rfx <- lme4::ranef(merMod, condVar=TRUE)
  out <- data.frame(groupFctr = NA, term = NA, estimate = NA,
                      std.error = NA, ER = NA, pctER = NA)

  if(!is.null(groupFctr)){
    groupFctr <- groupFctr
  } else{
    groupFctr <- names(rfx)
  }
  out <- data.frame(groupFctr = NA, groupLevel = NA, term = NA,
                    estimate = NA, std.error = NA,
                    ER = NA, pctER = NA)

  for(i in groupFctr){
      rfx.names <- rownames(rfx[[i]])
      n.grps <- length(rfx.names)
      n.terms <- length(rfx[[i]])
      if(!is.null(term)){
        termIdx <- term
      } else{
        termIdx <- names(rfx[[i]])
      }
      for(j in termIdx){
        if (grepl("[iI]nt[a-z]*", j) && is.na(match(j, names(rfx[[i]])))) {
           j <- "(Intercept)"
        }

        term.idx <- grep(j, names(rfx[[i]]), fixed=TRUE)
        theta <- rfx[[i]][,term.idx]
        var.theta <- attr(rfx[[i]], which="postVar")[term.idx, term.idx, 1:n.grps]
        #Calculate Expected Rank which is the sum of the probabilities that group i is greater than all
        #other groups j (assuming normal distribution of random effects)
        ER <- pctER <- rep(NA, n.grps)
        for (k in 1:n.grps) {
          ER[k] <- 1 + sum(pnorm((theta[k]-theta[-k]) / sqrt(var.theta[k] + var.theta[-k])))
        }
        #Calculated percentile expected rank ... the version of the formula I am using is
        #the percentage of groups that are ranked **equal to or less than** the selected
        #group ... if we just wanted percentage ranked less than then remove the 0.5
        pctER <- round(100 * (ER - 0.5) / n.grps)

        tmp <- data.frame(groupFctr = i, groupLevel = rfx.names, term = j,
                          estimate = theta, std.error = var.theta,
                          ER = ER, pctER = pctER)
        out <- rbind(out, tmp)
      }
  }

  out <- out[-1, ]
  # Avoid parentheses in parameter names
  out$term <- gsub("(Intercept)", "_Intercept", out$term, fixed = TRUE)
  return(out)
}
