#						AGREE.COEFF2.R
#					     (March 26, 2014)
# Description: This script file contains a series of R functions for computing 
# various agreement coefficients for 2 raters when the input data file is in 
# the form of 2x2 contingency table showing the distribution of subjects by 
# rater, and by category.
# 
# Author: Kilem L. Gwet, Ph.D.
#

#==============================================================================
# kappa2.table: Cohen's kappa (Cohen(1960)) coefficient and its standard error 
# for 2 raters when input dataset is a contingency table 
# -------------
# The input data "ratings" is a qxq contingency table showing the distribution 
# of subjects by rater, when q is the number of categories.
#==============================================================================

#' Cohen's kappa
#' 
#' Computes Cohen's kappa coefficient and standard error for two raters when 
#' data is a contingency table.
#' 
#' @param ratings a \code{q x q} contingency table showing the distribution of
#' subjects by rater where \code{q} is the number of categories
#' @param weights a vector of weights the same length as the number of 
#' categories
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @export

kappa2.table <- function(ratings, weights = diag(ncol(ratings)),
                         conflev = 0.95, N = Inf, print=TRUE) {
  if (dim(ratings)[1] != dim(ratings)[2])
  	stop('The contingency table should have the same number of rows and columns!') 
  n <- sum(ratings)               # number of subjects
  f <- n/N                        # final population correction  
  q <- ncol(ratings)              # number of categories 
  pa <- sum(weights * ratings/n)  # percent agreement
 
  pk. <- (ratings %*% rep(1, q)) / n
  p.l <- t((t(rep(1, q)) %*% ratings) / n)
  pe <- sum(weights*(pk. %*% t(p.l)))
  kappa <- (pa - pe) / (1 - pe)   # weighted kappa
  
	## 2 raters special case variance

  pkl <- ratings / n
  pb.k <- weights %*% p.l
  pbl. <- t(weights) %*% pk.
  sum1 <- 0
  for (k in 1:q) {
    for (l in 1:q) {
      sum1 <- sum1 + pkl[k,l] * (weights[k,l]-(1-kappa)*(pb.k[k] + pbl.[l]))^2
    }
  }
  var.kappa <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-kappa)*pe)^2)
  ## kappa standard error
  stderr <- sqrt(var.kappa)
  p.value <- 2*(1-pt(kappa/stderr,n-1))
  
  lcb <- kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if (print) {
    cat("Cohen's Kappa Coefficient\n")
    cat('=========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Kappa coefficient:',kappa,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 kappa = kappa,
                 stderr = stderr,
                 p.value = p.value))
}

# scott2.table: Scott's pi coefficient (Scott(1955)) and its standard error for 
# 2 raters when input dataset is a contingency table 
# -------------
# The input data "ratings" is a qxq contingency table showing the distribution
# of subjects by rater, when q is the number of categories.
#==============================================================================

#' Scott's pi coefficient
#' 
#' Computes Scott's pi coefficient and standard error for two raters when 
#' data is a contingency table.
#' 
#' @param ratings a \code{q x q} contingency table showing the distribution of
#' subjects by rater where \code{q} is the number of categories
#' @param weights a vector of weights the same length as the number of 
#' categories
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @export

scott2.table <- function(ratings, weights = diag(ncol(ratings)),
                         conflev = 0.95, N = Inf, print = TRUE) {
  if(dim(ratings)[1] != dim(ratings)[2])
  	stop('The contingency table should have the same number of rows and columns!') 
  
  n <- sum(ratings)              # number of subjects
  f <- n/N                       # final population correction  
  q <- ncol(ratings)             # number of categories 
  pa <- sum(weights * ratings/n) # percent agreement
 
  pk. <- (ratings %*% rep(1,q))/n
  p.l <- t((t(rep(1,q)) %*% ratings)/n)
  pi.k <- (pk.+p.l)/2
  pe <- sum(weights*(pi.k%*%t(pi.k)))
  scott <- (pa - pe)/(1 - pe)     # weighted scott's pi coefficint
  
	# 2 raters special case variance

  pkl <- ratings/n	         # p_{kl}	
  pb.k <- weights %*% p.l    # \ov{p}_{+k}
  pbl. <- t(weights) %*% pk. # \ov{p}_{l+}
  pbk  <- (pb.k + pbl.)/2    # \ov{p}_{k}
  sum1 <- 0
  for (k in 1:q) {
    for (l in 1:q){
      sum1 <- sum1 + pkl[k,l] * (weights[k,l]-(1-scott)*(pbk[k] + pbk[l]))^2
    }
  }
  
  var.scott <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-scott)*pe)^2)
  ## Scott's standard error
  stderr <- sqrt(var.scott)
  p.value <- 2*(1-pt(scott/stderr,n-1))
  
  lcb <- scott - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,scott + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print) {
    cat("Scott's Pi Coefficient\n")
    cat('======================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Scott coefficient:',scott,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 scott = scott,
                 stderr = stderr,
                 p.value = p.value))
}

# gwet.ac1.table: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard 
# error for 2 raters when input dataset is a contingency table 
# -------------
# The input data "ratings" is a qxq contingency table showing the distribution 
# of subjects by rater, when q is the number of categories.
#==============================================================================

#' Gwet's AC1/AC2 coefficient
#' 
#' Computes Gwet's AC1/AC2 coefficient and standard error for two raters when 
#' data is a contingency table.
#' 
#' @param ratings a \code{q x q} contingency table showing the distribution of
#' subjects by rater where \code{q} is the number of categories
#' @param weights a vector of weights the same length as the number of 
#' categories
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @export

gwet.ac1.table <- function(ratings, weights = diag(ncol(ratings)), 
                           conflev = 0.95, N = Inf, print = TRUE) {
  if(dim(ratings)[1] != dim(ratings)[2])
  	stop('The contingency table should have the same number of rows and columns!')
  
  n <- sum(ratings)              # number of subjects
  f <- n/N                       # final population correction  
  q <- ncol(ratings)             # number of categories 
  pa <- sum(weights * ratings/n) # percent agreement
 
  pk. <- (ratings %*% rep(1,q))/n
  p.l <- t((t(rep(1,q)) %*% ratings)/n)
  pi.k <- (pk. + p.l)/2
  tw <- sum(weights)
  pe <- tw * sum(pi.k * (1-pi.k)) / (q * (q-1))
  # gwet's ac1/ac2 coefficint
  gwet.ac1 <- (pa - pe)/(1 - pe) 
  
	# calculation of variance - standard error - confidence interval - p-value

  pkl <- ratings/n	     #p_{kl}	
  sum1 <- 0
  for (k in 1:q) {
    for (l in 1:q) {
      sum1 <- sum1 + pkl[k,l] * (weights[k,l]-2*(1-gwet.ac1)*tw*
                                   (1-(pi.k[k] + pi.k[l])/2)/(q*(q-1)))^2
    }
  }
  
  var.gwet <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-gwet.ac1)*pe)^2)
  # ac1's standard error
  stderr <- sqrt(var.gwet)
  p.value <- 2*(1-pt(gwet.ac1/stderr,n-1))
  
  # confidence bounds
  lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1) 
  ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Gwet's AC1/AC2 Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('AC1/AC2 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 gwet.ac1 = gwet.ac1,
                 stderr = stderr,
                 p.value = p.value))
}

# bp2.table: Brennan-Prediger coefficient (Brennan & Prediger (1981)) and its
# standard error for 2 raters when input dataset is a contingency table 
# -------------
# The input data "ratings" is a qxq contingency table showing the distribution
# of subjects by rater, when q is the number of categories.
#==============================================================================

#' Brennan-Prediger coefficient
#' 
#' Computes Brennan-Prediger coefficient and standard error for two raters when
#' data is a contingency table.
#' 
#' @param ratings a \code{q x q} contingency table showing the distribution of
#' subjects by rater where \code{q} is the number of categories
#' @param weights a vector of weights the same length as the number of 
#' categories
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @export

bp2.table <- function(ratings, weights = diag(ncol(ratings)),
                      conflev = 0.95, N = Inf, print = TRUE) {
  if (dim(ratings)[1] != dim(ratings)[2])
  	stop('The contingency table should have the same number of rows and columns!')
  
  n <- sum(ratings)              # number of subjects
  f <- n/N                       # final population correction  
  q <- ncol(ratings)             # number of categories 
  pa <- sum(weights * ratings/n) # percent agreement
 
  tw <- sum(weights)
  pe <- tw/(q^2)
  # Brennan-Prediger coefficint
  bp.coeff <- (pa - pe)/(1 - pe)
  
	# calculation of variance - standard error - confidence interval - p-value

  pkl <- ratings/n	     #p_{kl}	
  sum1 <- 0
  for (k in 1:q) {
    for (l in 1:q) {
      sum1 <- sum1 + pkl[k,l] * weights[k,l]^2
    }
  }
  
  var.bp <- ((1-f)/(n*(1-pe)^2)) * (sum1 - pa^2)
  # bp's standard error
  stderr <- sqrt(var.bp)
  p.value <- 2*(1-pt(bp.coeff/stderr,n-1))
  
  ## confidence bounds
  lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Brennan-Prediger Coefficient\n")
    cat('============================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('B-P coefficient:',bp.coeff,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 bp.coeff = bp.coeff,
                 stderr = stderr,
                 p.value = p.value))
}

# krippen2.table: Krippen's alpha coefficient (Scott(1955)) and its standard 
# error for 2 raters when input dataset is a contingency table 
# -------------
# The input data "ratings" is a qxq contingency table showing the distribution 
# of subjects by rater, when q is the number of categories.
#==============================================================================

#' Krippendorff's alpha coefficient
#' 
#' Computes Krippendorff's alpha coefficient and standard error for two raters 
#' when data is a contingency table.
#' 
#' @param ratings a \code{q x q} contingency table showing the distribution of
#' subjects by rater where \code{q} is the number of categories
#' @param weights a vector of weights the same length as the number of 
#' categories
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @export

krippen2.table <- function(ratings, weights = diag(ncol(ratings)),
                           conflev = 0.95, N = Inf, print = TRUE) {
  if (dim(ratings)[1] != dim(ratings)[2])
  	stop('The contingency table should have the same number of rows and columns!') 
  
  n <- sum(ratings)  # number of subjects
  f <- n/N           # final population correction  
  q <- ncol(ratings) # number of categories 
  epsi = 1/(2*n)
  pa0 <- sum(weights * ratings/n)
  # percent agreement
  pa <- (1-epsi)*pa0 + epsi 
 
  pk. <- (ratings %*% rep(1,q))/n
  p.l <- t((t(rep(1,q)) %*% ratings)/n)
  pi.k <- (pk.+p.l)/2
  pe <- sum(weights*(pi.k%*%t(pi.k)))
  # weighted Krippen's alpha coefficint
  kripp.coeff <- (pa - pe)/(1 - pe) 
  
	# calculating variance

  pkl <- ratings/n	     #p_{kl}	
  pb.k <- weights %*% p.l    #\ov{p}_{+k}
  pbl. <- t(weights) %*% pk. #\ov{p}_{l+}
  pbk  <- (pb.k + pbl.)/2    #\ov{p}_{k}
  sum1 <- 0
  for (k in 1:q) {
    for (l in 1:q) {
      sum1 <- sum1 + pkl[k,l] * ((1-epsi)*weights[k,l]-(1-kripp.coeff)*
                                   (pbk[k] + pbk[l]))^2
    }
  }
  
  var.kripp <- ((1-f)/(n*(1-pe)^2)) * (sum1 - ((1-epsi)*
                                                 pa0-2*(1-kripp.coeff)*pe)^2)
  # Kripp. alpha's standard error
  stderr <- sqrt(var.kripp)
  p.value <- 2*(1-pt(kripp.coeff/stderr,n-1))
  
  ## confidence bounds
  lcb <- kripp.coeff - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,kripp.coeff + stderr*qt(1-(1-conflev)/2,n-1))
  if(print) {
    cat("Krippendorff's alpha Coefficient\n")
    cat('================================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Alpha coefficient:',kripp.coeff,'Standard error:',stderr,'\n')
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 kripp.coeff = kripp.coeff,
                 stderr = stderr,
                 p.value = p.value))
}