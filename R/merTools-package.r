#' merTools: Provides methods for extracting and exploring results from merMod
#' objects in the lme4 package.
#'
#' The merTools package contains convenience tools for extracting useful
#' information from and exploring the implications of merMod objects created by
#' the lme4 package.  These convenience functions are especially useful for
#' merMod objects that take a long time to estimate due to their complexity or
#' because they are estimated on very large samples.
#'
#' See the vignettes for usage examples
#'
#' @section merMod extraction/utility functions:
#'
#' \itemize{
#'   \item \code{\link{fastdisp}}
#'   \item \code{\link{superFactor}}
#'   \item \code{\link{REextract}}
#'   \item \code{\link{REsim}}
#'   \item \code{\link{FEsim}}
#'   \item \code{\link{RMSE.merMod}}
#'   \item \code{\link{thetaExtract}}
#'   \item \code{\link{REquantile}}
#' }
#'
#' @section merMod exploration functions:
#'
#' \itemize{
#'   \item \code{\link{plotREsim}}
#'   \item \code{\link{plotFEsim}}
#'   \item \code{\link{draw}}
#'   \item \code{\link{wiggle}}
#'   \item \code{\link{subBoot}}
#'   \item \code{\link{predictInterval}}
#'   \item \code{\link{expectedRank}}
#'   \item \code{\link{REimpact}}
#'   \item \code{\link{shinyMer}}
#' }
#'
#' @name merTools
#' @docType package
#' @aliases merTools merTools-package
NULL

#' A subset of data from the 1982 High School and Beyond survey used as examples for HLM software
#' @description A key example dataset used for examples in the HLM software manual.
#' Included here for use in replicating HLM analyses in R.
#' @format A data frame with 7,185 observations on the following 8 variables.
#' \describe{
#'  \item{\code{schid}}{a numeric vector, 160 unique values}
#'  \item{\code{mathach}}{a numeric vector for the performance on a standardized math assessment}
#'  \item{\code{female}}{a numeric vector coded 0 for male and 1 for female}
#'  \item{\code{ses}}{a numeric measure of student socio-economic status}
#'  \item{\code{minority}}{a numeric vector coded 0 for white and 1 for non-white students}
#'  \item{\code{schtype}}{a numeric vector coded 0 for public and 1 for private schools}
#'  \item{\code{meanses}}{a numeric, the average SES for each school in the data set}
#'  \item{\code{size}}{a numeric for the number of students in the school}
#' }
#' @details The data file used for this presentation is a subsample from the
#' 1982 High School and Beyond Survey and is used extensively in
#' Hierarchical Linear Models by Raudenbush and Bryk. It consists of 7,185 students
#' nested in 160 schools.
#' @source Data made available by UCLA Institute for Digital Research and Education
#' (IDRE) online: \url{http://www.ats.ucla.edu/stat/hlm/seminars/hlm6/mlm_hlm6_seminar.htm}
#' @references Stephen W. Raudenbush and Anthony S. Bryk (2002). Hierarchical
#' Linear Models: Applications and Data Analysis Methods (2nd ed.). SAGE.
#' @examples
#' data(hsb)
#' head(hsb)
"hsb"


