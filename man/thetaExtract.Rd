% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subBoot.R
\name{thetaExtract}
\alias{thetaExtract}
\title{Extract theta parameters from a merMod model}
\usage{
thetaExtract(merMod)
}
\arguments{
\item{merMod}{a valid merMod object}
}
\value{
a vector of the covariance, theta, parameters from a \code{\link{merMod}}
}
\description{
A convenience function that returns the theta parameters for a
\code{\link{merMod}} object.
}
\examples{
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
thetaExtract(fm1) #(a numeric vector of the covariance parameters)
}
\seealso{
merMod
}
