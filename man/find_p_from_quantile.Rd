% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility-functions.R
\name{find_p_from_quantile}
\alias{find_p_from_quantile}
\title{Find logarithmic distribution parameter from a quantile value}
\usage{
find_p_from_quantile(q, conf = 0.95)
}
\arguments{
\item{q}{An integer equal to, or greater than 1. Quantile value.}

\item{conf}{A probability. Cumulative probability of the
logarithmic distribution at \code{q}.}
}
\value{
parameter for the logarithmic distribution
}
\description{
Find logarithmic distribution parameter from a quantile value
}
