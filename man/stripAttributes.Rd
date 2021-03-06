% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merData.R
\name{stripAttributes}
\alias{stripAttributes}
\title{Remove attributes from a data.frame}
\usage{
stripAttributes(data)
}
\arguments{
\item{data}{a data.frame}
}
\value{
a data frame with variable names cleaned to remove all attributes except for
names, row.names, and class
}
\description{
Strips attributes off of a data frame that come with a merMod model.frame
}
