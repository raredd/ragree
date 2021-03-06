% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwet_agree.coeff3.dist.r
\name{bp.coeff.dist}
\alias{bp.coeff.dist}
\title{Brennan-Prediger coefficient}
\usage{
bp.coeff.dist(ratings, weights = "unweighted", conflev = 0.95,
  N = Inf, print = TRUE)
}
\arguments{
\item{ratings}{an \code{n x q} matrix showing the number of raters by
subject and category where \code{n} is the number of subjects and \code{q}
is the number of categories}

\item{weights}{optional weighting}

\item{conflev}{confidence level}

\item{N}{used as denominator in finite population correction}

\item{print}{logical; if \code{TRUE}, prints a summary of the agreement}
}
\description{
Computes Brennan-Prediger coefficient and standard error for multiple 
raters when data is an \code{n x q} matrix representing the distribution of
raters be subject and by category.
}
\details{
A typical entry associated with a subject and a category, represents the 
number of raters who classified the subject into the specified category. 
Excludes all subjects that are not rated by any rater.
}
\references{
Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, 
misuses, and alternatives." Educational and Psychological Measurement, 41, 
687-699.
}
\author{
Kilem L. Gwet
}
