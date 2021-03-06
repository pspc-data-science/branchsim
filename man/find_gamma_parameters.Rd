% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMJ-functions.R
\name{find_gamma_parameters}
\alias{find_gamma_parameters}
\title{Get the shape and rate parameter of a gamma distribution given an integration interval and mean.}
\usage{
find_gamma_parameters(m = 5.5, a = 0, b = 11.5, int_value = 0.975)
}
\arguments{
\item{m}{The mean of the gamma distribution.}

\item{a}{The lower limit of integration. Default is a=0}

\item{b}{The upper limit of intefration. Default is b=11.5}

\item{int_value}{The confidence level \eqn{\int_a^b f(x)\,dx = int_value}. Default is int_value = .975}

\item{P}{The parameter of the logarithmic distribution for the number of infected during an event.
Default p=0.5}
}
\value{
The shape and rate parameters of the gamma distribution.
}
\description{
Get the shape and rate parameter of a gamma distribution given an integration interval and mean.
}
