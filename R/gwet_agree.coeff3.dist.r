#					AGREE.COEFF3.DIST.R
#				 	 (March 26, 2014)
# Description: This script file contains a series of R functions for computing
# various agreement coefficients for multiple raters (2 or more) when the input
# data file is in the form of nxq matrix or data frame showing the count of 
# raters by subject and by category. 
#
# That is n = number of subjects, and q = number of categories.
#
# A typical table entry (i,k) represents the number of raters who classified 
# subject i into category k. 
# 
# Author: Kilem L. Gwet, Ph.D.

#==============================================================================
# gwet.ac1.dist: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error 
# for multiple raters when input dataset is a nxq matrix representing the 
# distribution of raters by subject and by category. 
#-------------
# The input data "ratings" is an nxq matrix showing the number of raters by 
# subject and category. A typical entry associated with a subject and a 
# category, represents the number of raters who classified the subject into the
# specified category. Exclude all subjects that are not rated by any rater.
#
# Bibliography:
# Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in 
#     the presence of high agreement." British Journal of Mathematical and 
#     Statistical Psychology, 61, 29-48.
#==============================================================================

#' Gwet's AC1/AC2 coefficient
#' 
#' Computes Gwet's AC1/AC2 coefficient and standard error for multiple raters
#' when data is an \code{n x q} matrix representing the distribution of 
#' raters be subject and by category.
#' 
#' A typical entry associated with a subject and a category, represents the 
#' number of raters who classified the subject into the specified category. 
#' Excludes all subjects that are not rated by any rater.
#' 
#' @param ratings an \code{n x q} matrix showing the number of raters by
#' subject and category where \code{n} is the number of subjects and \code{q}
#' is the number of categories
#' @param weights optional weighting
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references Gwet, K. L. (2008). "Computing inter-rater reliability and its
#' variance in the presence of high agreement." British Journal of Mathematical
#' and Statistical Psychology, 61, 29 - 48.
#' @export

gwet.ac1.dist <- function(ratings, weights = "unweighted",
                          conflev = 0.95, N = Inf, print = TRUE) { 
  agree.mat <- as.matrix(ratings) 
  n <- nrow(agree.mat) # number of subjects
  q <- ncol(agree.mat) # number of categories
  f <- n/N             # final population correction 

  # creating the weights matrix

  if (is.character(weights)) {
     weights.mat <- diag(q)
  } else weights.mat = as.matrix(weights)
  
  agree.mat.w <- t(weights.mat %*% t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat %*% rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1)) %*% rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n)) %*% (agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  gwet.ac1 <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next ratio 
  # calculation always possible.
  den.ivec <- den.ivec - (den.ivec==0)
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-gwet.ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - gwet.ac1)^2)
  stderr <- sqrt(var.ac1)# ac1's standard error
  p.value <- 2*(1-pt(gwet.ac1/stderr,n-1))
  
  # confidence bounds
  lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Gwet's AC1/AC2 Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    if (!is.character(weights)) {
	  cat('AC2 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
	  cat('Weights:\n')
	  write.table(weights,row.names=FALSE,col.names=FALSE)
        cat('\n')
    } else
      cat('AC1 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 gwet.ac1 = gwet.ac1,
                 stderr = stderr,
                 p.value = p.value))
}

#==============================================================================
# fleiss.kappa.dist: This function computes Fleiss' generalized kappa 
# coefficient (see Fleiss(1971)) and its standard error for 3 raters or more 
# when input dataset is a nxq matrix representing the distribution of raters by
# subject and by category. 
# -------------
# The input data "ratings" is an nxq matrix showing the number of raters by 
# subject and category. A typical entry associated with a subject and a 
# category, represents the number of raters who classified the subject into the
# specified category. Exclude all subjects that are not rated by any rater.
#
# Bibliography:
# Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John 
#     Wiley & Sons.
#==============================================================================

#' Fleiss's generalized kappa coefficient
#' 
#' Computes Fleiss's kappa coefficient and standard error for multiple raters
#' when data is an \code{n x q} matrix representing the distribution of 
#' raters be subject and by category.
#' 
#' A typical entry associated with a subject and a category, represents the 
#' number of raters who classified the subject into the specified category. 
#' Excludes all subjects that are not rated by any rater.
#' 
#' @param ratings an \code{n x q} matrix showing the number of raters by
#' subject and category where \code{n} is the number of subjects and \code{q}
#' is the number of categories
#' @param weights optional weighting
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references Fleiss, J. L. (1981). Statistical Methods for Rates and 
#' Proportions. John Wiley & Sons.
#' @export

fleiss.kappa.dist <- function(ratings, weights = "unweighted",
                              conflev = 0.95, N = Inf, print = TRUE) { 
  agree.mat <- as.matrix(ratings) 
  n <- nrow(agree.mat) # number of subjects
  q <- ncol(agree.mat) # number of categories
  f <- n/N             # final population correction 

  # creating the weights matrix

  if (is.character(weights)) {
     weights.mat<-diag(q)
  } else weights.mat= as.matrix(weights)
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating fleiss's generalized kappa coefficient

  ri.vec <- agree.mat %*% rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1)) %*% rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n)) %*% (agree.mat/(ri.vec %*% t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec %*% t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next 
  # ratio calculation always possible.
  den.ivec <- den.ivec - (den.ivec==0)
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat %*% pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.fleiss <- ((1-f)/(n*(n-1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
  stderr <- sqrt(var.fleiss)# kappa's standard error
  p.value <- 2*(1-pt(fleiss.kappa/stderr,n-1))
  
  # confidence bounds
  lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Fleiss' Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Fleiss kappa coefficient:',fleiss.kappa,'Standard error:',stderr,'\n')
    if (!is.character(weights)) {
	  cat('Weights:\n')
	  write.table(weights,row.names=FALSE,col.names=FALSE)
    cat('\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(c(pa,pe,fleiss.kappa,stderr,p.value))
  invisible(list(pa = pa,
                 pe = pe,
                 fleiss.kappa = fleiss.kappa,
                 stderr = stderr,
                 p.value = p.value))
}

#==============================================================================
# krippen.alpha.dist: This function computes Krippendorff's alpha coefficient 
# (see Krippendorff(1970, 1980)) and its standard error for 3 raters or more 
# when input dataset is a nxq matrix representing the distribution of raters by
# subject and by category. 
# 
# The input data "ratings" is an nxq matrix showing the number of raters by 
# subject and category. A typical entry associated with a subject and a 
# category, represents the number of raters who classified the subject into the
# specified category. Exclude all subjects that are not rated by any rater.
#-------------
# The algorithm used to compute krippendorff's alpha is very different from 
# anything that was published on this topic. Instead, it follows the equations 
# presented by K. Gwet (2010).
#
# Bibliography:
# Gwet, K. (2012). Handbook of Inter-Rater Reliability: The Definitive Guide to
# Measuring the Extent of Agreement Among Multiple Raters, 3rd Edition. 
# Advanced Analytics, LLC; 3rd edition (March 2, 2012)
# Krippendorff (1970). "Bivariate agreement coefficients for reliability of 
# data." Sociological Methodology,2,139-150
# Krippendorff (1980). Content analysis: An introduction to its methodology 
# (2nd ed.), New-bury Park, CA: Sage.
#==============================================================================

#' Krippendorff's alpha coefficient
#' 
#' Computes Krippendorff's alpha coefficient and standard error for multiple 
#' raters when data is an \code{n x q} matrix representing the distribution of
#' raters be subject and by category.
#' 
#' A typical entry associated with a subject and a category, represents the 
#' number of raters who classified the subject into the specified category. 
#' Excludes all subjects that are not rated by any rater.
#' 
#' The algorithm used to compute Krippendorff's alpha is very different from 
#' anything that was published on this topic. Instead, it follows the equations
#' presented by K. Gwet (2010).
#' 
#' @param ratings an \code{n x q} matrix showing the number of raters by
#' subject and category where \code{n} is the number of subjects and \code{q}
#' is the number of categories
#' @param weights optional weighting
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references
#' Gwet, K. (2012). Handbook of Inter-Rater Reliability: the Definitive Guide 
#' to Measuring the Extent of Agreement among Multiple Raters, 3rd Edition.
#' Advanced Analytics, LLC; 3rd edition (March 2, 2012).
#' 
#' Krippendorff (1970). "Bivariate agreement coefficients for reliability of 
#' data." Sociological Methodology, 2, 139-150.
#' 
#' Krippendorff (1980). Content analysis: An introduction to its methodology 
#' (2nd ed.), New-bury Park, CA: Sage.
#' @export

krippen.alpha.dist <- function(ratings, weights = "unweighted",
                               conflev = 0.95, N = Inf, print = TRUE) { 
  agree.mat <- as.matrix(ratings) 
  n <- nrow(agree.mat) # number of subjects
  q <- ncol(agree.mat) # number of categories
  f <- n/N             # final population correction 

  # creating the weights matrix
  if (is.character(weights)) {
     weights.mat <- diag(q)
  } else weights.mat= as.matrix(weights)
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating krippendorff's alpha coefficient
  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat<-agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  pa <- (1-epsi)* sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.mean*(ri.vec-1)
  pa.ivec <- sum.q/den.ivec
  pa.v <- mean(pa.ivec)
  pa.ivec <- (1-epsi)*(pa.ivec-pa.v*(ri.vec-ri.mean)/ri.mean) + epsi

  krippen.ivec <- (pa.ivec-pe)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec

  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.mean - sum(pi.vec) * 
    (ri.vec-ri.mean)/ri.mean
  krippen.ivec.x <- krippen.ivec - (1-krippen.alpha) * (pe.ivec-pe)/(1-pe)
  
  var.krippen <- ((1-f)/(n*(n-1))) * sum((krippen.ivec.x - krippen.alpha)^2)
  # alpha's standard error
  stderr <- sqrt(var.krippen)
  p.value <- 2*(1-pt(krippen.alpha/stderr,n-1))
  
  # confidence bounds
  lcb <- krippen.alpha - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,krippen.alpha + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Krippendorff's Alpha Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Krippendorff alpha coefficient:',krippen.alpha,'Standard error:',
        stderr,'\n')
    if (!is.character(weights)) {
	  cat('Weights:\n')
	  write.table(weights,row.names=FALSE,col.names=FALSE)
    cat('\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 krippen.alpha = krippen.alpha,
                 stderr = stderr,
                 p.value = p.value))
}

#==============================================================================
#bp.coeff.dist: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and
# its standard error for multiple raters when input dataset is a nxq matrix 
# representing the distribution of raters by subject and by category. 
# 
# The input data "ratings" is an nxq matrix showing the number of raters by 
# subject and category. A typical entry associated with a subject and a 
# category, represents the number of raters who classified the subject into the
# specified category. Exclude all subjects that are not rated by any rater.
#--------------------------------------------
# Bibliography:
# Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, 
#     misuses, and alternatives." Educational and Psychological Measurement, 41, 
#     687-699.
#==============================================================================

#' Brennan-Prediger coefficient
#' 
#' Computes Brennan-Prediger coefficient and standard error for multiple 
#' raters when data is an \code{n x q} matrix representing the distribution of
#' raters be subject and by category.
#' 
#' A typical entry associated with a subject and a category, represents the 
#' number of raters who classified the subject into the specified category. 
#' Excludes all subjects that are not rated by any rater.
#' 
#' @param ratings an \code{n x q} matrix showing the number of raters by
#' subject and category where \code{n} is the number of subjects and \code{q}
#' is the number of categories
#' @param weights optional weighting
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references
#' Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, 
#' misuses, and alternatives." Educational and Psychological Measurement, 41, 
#' 687-699.
#' @export

bp.coeff.dist <- function(ratings, weights = "unweighted", 
                          conflev = 0.95, N = Inf, print = TRUE) {
  agree.mat <- as.matrix(ratings) 
  n <- nrow(agree.mat) # number of subjects
  q <- ncol(agree.mat) # number of categories
  f <- n/N             # final population correction 

  # creating the weights matrix
  if (is.character(weights)) {
     weights.mat<-diag(q)
  } else weights.mat= as.matrix(weights)
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next 
  # ratio calculation always possible.
  den.ivec <- den.ivec - (den.ivec==0)
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  bp.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  var.bp <- ((1-f)/(n*(n-1))) * sum((bp.ivec - bp.coeff)^2)
  # BP's standard error
  stderr <- sqrt(var.bp)
  p.value <- 2*(1-pt(bp.coeff/stderr,n-1))
  
  # confidence bounds
  lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1))
  if (print) {
    cat("Brennan-Prediger Coefficient\n")
    cat('============================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('B-P coefficient:',bp.coeff,'Standard error:',stderr,'\n')
    if (!is.character(weights)) {
	  cat('Weights:\n')
	  write.table(weights,row.names=FALSE,col.names=FALSE)
    cat('\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 bp.coeff = bp.coeff,
                 stderr = stderr,
                 p.value = p.value))
}
