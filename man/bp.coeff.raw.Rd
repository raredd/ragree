% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwet_agree.coeff3.raw.r
\name{bp.coeff.raw}
\alias{bp.coeff.raw}
\title{Brennan-Prediger kappa coefficient}
\usage{
bp.coeff.raw(ratings, weights = "unweighted", conflev = 0.95,
  N = Inf, print = TRUE)
}
\arguments{
\item{ratings}{an \code{n x r} matrix of raw alphanumeric ratings from
\code{n} subjects and \code{r} raters}

\item{weights}{optional weighting; see details}

\item{conflev}{confidence level}

\item{N}{used as denominator in finite population correction}

\item{print}{logical; if \code{TRUE}, prints a summary of the agreement}
}
\description{
Computes Brennan-Prediger kappa coefficient and standard error for multiple
raters when data is an \code{n x r} matrix of alphanumeric ratings from 
\code{n} subjects and \code{r} raters, excluding all subjects that are not 
rated by any rater.
}
\details{
\code{weight} is an option matrix of weights or one of \code{"quadratic"},
\code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
\code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
matrix.
}
\references{
Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, 
misuses, and alternatives." Educational and Psychological Measurement, 41,
687-699.
}
\author{
Kilem L. Gwet
}
