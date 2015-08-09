#								AGREE.COEFF3.RAW.R
#						 		 (March 26, 2014)
# Description: This script file contains a series of R functions for computing
# various agreement coefficients for multiple raters (2 or more) when the 
# input data file is in the form of nxr matrix or data frame showing the actual
# ratings each rater (column) assigned to each subject (in row). 
# 
# That is n = number of subjects, and r = number of raters.
# 
# A typical table entry (i,g) represents the rating associated with subject i 
# and rater g. 
#
# Author: Kilem L. Gwet, Ph.D.

trim <- function(x) gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)

#==============================================================================
# gwet.ac1.raw: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error 
# for multiple raters when input dataset is a nxr matrix of alphanumeric 
# ratings from n subjects and r raters 
# -------------
# The input data "ratings" is a nxr data frame of raw alphanumeric ratings from
# n subjects and r raters. Exclude all subjects that are not rated by any rater
# 
# Bibliography:
# Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in 
#     the presence of high agreement." British Journal of Mathematical and 
#     Statistical Psychology, 61, 29-48.
#==============================================================================

#' Gwet's AC1/AC2 coefficient
#' 
#' Computes Gwet's AC1/AC2 coefficient and standard error for multiple raters
#' when data is an \code{n x r} matrix of alphanumeric ratings from \code{n}
#' subjects and \code{r} raters, excluding all subjects that are not rated by 
#' any rater.
#' 
#' \code{weight} is an option matrix of weights or one of \code{"quadratic"},
#' \code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
#' \code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
#' matrix.
#' 
#' @param ratings an \code{n x r} matrix of raw alphanumeric ratings from
#' \code{n} subjects and \code{r} raters
#' @param weights optional weighting; see details
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references Gwet, K. L. (2008). "Computing inter-rater reliability and its
#' variance in the presence of high agreement." British Journal of Mathematical
#' and Statistical Psychology, 61, 29 - 48.
#' @export

gwet.ac1.raw <- function(ratings, weights = "unweighted", 
                         conflev = 0.95, N = Inf, print = TRUE) { 
  
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N               # final population correction 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else {
    #trim vector elements to remove leading and trailing blanks
    categ.init <- trim(categ.init)
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  } else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution 
  # of raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for (k in 1:q) {
    if (is.numeric(ratings.mat)) {
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    } else
      agree.mat[,k] <- (ratings.mat==categ[k]) %*% rep(1,r)
  }
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
  if(print) {
    cat("Gwet's AC1/AC2 Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    if (weights=="unweighted")
      cat('AC1 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
    else {
      cat('AC2 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
      if (!is.numeric(weights)) {
        cat('Weights: ', weights,'\n')
      } else
        cat('Weights: Custom Weights\n')
    }
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
# fleiss.kappa.raw: This function computes Fleiss generalized kappa coefficient
# (see Fleiss(1971)) and its standard error for 3 raters or more when input 
# dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters.
#-------------
# The input data "ratings" is a nxr data frame of raw alphanumeric ratings
# from n subjects and r raters. Exclude all subjects that are not rated by any 
# rater.
# 
# Bibliography:
# Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John 
# Wiley & Sons.
#==============================================================================

#' Fleiss's generalized kappa coefficient
#' 
#' Computes Fleiss's kappa coefficient and standard error for multiple raters
#' when data is an \code{n x r} matrix of alphanumeric ratings from \code{n}
#' subjects and \code{r} raters, excluding all subjects that are not rated by 
#' any rater.
#' 
#' \code{weight} is an option matrix of weights or one of \code{"quadratic"},
#' \code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
#' \code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
#' matrix.
#' 
#' @param ratings an \code{n x r} matrix of raw alphanumeric ratings from
#' \code{n} subjects and \code{r} raters
#' @param weights optional weighting; see details
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references 
#' Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John 
#' Wiley and Sons.
#' @export

fleiss.kappa.raw <- function(ratings, weights = "unweighted", 
                             conflev = 0.95, N = Inf, print = TRUE) {
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N               # final population correction 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
     categ <- sort(as.vector(na.omit(categ.init)))
  } else {
    # trim vector elements to remove leading and trailing blanks
    categ.init <- trim(categ.init)
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of 
  # raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for (k in 1:q) {
    if (is.numeric(ratings.mat)) {
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k %*% rep(1,r) 
      } else
        agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat %*% t(agree.mat))

  # calculating fleiss's generalized kappa coefficient
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n)) %*% (agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec %*% t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next ratio 
  # calculation always possible.
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
  # kappa's standard error
  stderr <- sqrt(var.fleiss)
  p.value <- 2*(1-pt(fleiss.kappa/stderr,n-1))
  
  # confidences bounds
  lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Fleiss' Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Fleiss kappa coefficient:',fleiss.kappa,'Standard error:',stderr,'\n')
    
    if (weights!="unweighted") {
      if (!is.numeric(weights)) {
        cat('Weights: ', weights,'\n')
      } else
        cat('Weights: Custom Weights\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(list(pa = pa,
                 pe = pe,
                 fleiss.kappa = fleiss.kappa,
                 stderr = stderr,
                 p.value = p.value))
}

#==============================================================================
# krippen.alpha.raw: This function computes Krippendorff's alpha coefficient 
# (see Krippendorff(1970, 1980)) and its standard error for 3 raters or more 
# when input dataset is a nxr matrix of alphanumeric ratings from n subjects 
# and r raters.
#-------------
# The algorithm used to compute krippendorff's alpha is very different from 
# anything that was published on this topic. Instead, it follows the equations
# presented by K. Gwet (2010)
# 
# The input data "ratings" is a nxr data frame of raw alphanumeric ratings
# from n subjects and r raters. Exclude all subjects that are not rated by 
# any rater.
# 
# Bibliography:
# Gwet, K. (2012). Handbook of Inter-Rater Reliability: The Definitive Guide 
#     to Measuring the Extent of Agreement Among Multiple Raters, 3rd Edition.
#     Advanced Analytics, LLC; 3rd edition (March 2, 2012).
# Krippendorff (1970). "Bivariate agreement coefficients for reliability of 
# data." Sociological Methodology,2,139-150.
# Krippendorff (1980). Content analysis: An introduction to its methodology 
# (2nd ed.), New-bury Park, CA: Sage.
#==============================================================================


#' Krippendorff's alpha coefficient
#' 
#' Computes Krippendorff's alpha coefficient and standard error for multiple 
#' raters when data is an \code{n x r} matrix of alphanumeric ratings from 
#' \code{n} subjects and \code{r} raters, excluding all subjects that are not 
#' rated by any rater.
#' 
#' \code{weight} is an option matrix of weights or one of \code{"quadratic"},
#' \code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
#' \code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
#' matrix.
#' 
#' The algorithm used to compute Krippendorff's alpha is very different from 
#' anything that was published on this topic. Instead, it follows the equations
#' presented by K. Gwet (2010).
#' 
#' @param ratings an \code{n x r} matrix of raw alphanumeric ratings from
#' \code{n} subjects and \code{r} raters
#' @param weights optional weighting; see details
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

krippen.alpha.raw <- function(ratings, weights = "unweighted", 
                              conflev = 0.95, N = Inf, print = TRUE) { 
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N               # final population correction 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else {
    #trim vector elements to remove leading and trailing blanks
    categ.init <- trim(categ.init)
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of 
  # raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)) {
    k.mis <-(ratings.mat==categ[k])
    in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    } else
      agree.mat[,k] <- (ratings.mat==categ[k]) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat %*% t(agree.mat))

  # calculating krippendorff's alpha coefficient
  ri.vec <- agree.mat %*% rep(1,q)
  agree.mat<-agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1)) %*% rep(1,q)
  pa <- (1-epsi)* sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi

  pi.vec <- t(t(rep(1/n,n)) %*% (agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.mean*(ri.vec-1)
  pa.ivec <- sum.q/den.ivec
  pa.v <- mean(pa.ivec)
  pa.ivec <- (1-epsi)*(pa.ivec-pa.v*(ri.vec-ri.mean)/ri.mean) + epsi

  krippen.ivec <- (pa.ivec-pe)/(1-pe)
  pi.vec.wk. <- weights.mat %*% pi.vec
  pi.vec.w.k <- t(weights.mat) %*% pi.vec

  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat %*% pi.vec.w)/ri.mean - sum(pi.vec) * 
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
    if (weights!="unweighted") {
      if (!is.numeric(weights)) {
        cat('Weights: ', weights,'\n')
	 } else
     cat('Weights: Custom Weights\n')
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
# conger.kappa.raw: Conger's kappa coefficient (see Conger(1980)) and its 
# standard error for multiple raters when input dataset is a nxr matrix of 
# alphanumeric ratings from n subjects and r raters 
# -------------
# The input data "ratings" is a nxr data frame of raw alphanumeric ratings
# from n subjects and r raters. Exclude all subjects that are not rated by any 
# rater.
#
# Bibliography:
# Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple 
# Raters," Psychological Bulletin, 88, 322-328.
#==============================================================================

#' Conger's kappa coefficient
#' 
#' Computes Conger's kappa coefficient and standard error for multiple 
#' raters when data is an \code{n x r} matrix of alphanumeric ratings from 
#' \code{n} subjects and \code{r} raters, excluding all subjects that are not 
#' rated by any rater.
#' 
#' \code{weight} is an option matrix of weights or one of \code{"quadratic"},
#' \code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
#' \code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
#' matrix.
#' 
#' @param ratings an \code{n x r} matrix of raw alphanumeric ratings from
#' \code{n} subjects and \code{r} raters
#' @param weights optional weighting; see details
#' @param conflev confidence level
#' @param N used as denominator in finite population correction
#' @param print logical; if \code{TRUE}, prints a summary of the agreement
#' 
#' @author Kilem L. Gwet
#' @references
#' Conger, A. J. (1980), "Integration and Generalization of Kappas for Multiple
#' Raters," Psychological Bulletin, 88, 322-328.
#' @export

conger.kappa.raw <- function(ratings, weights = "unweighted",
                             conflev = 0.95, N = Inf, print = TRUE) {
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N               # final population correction 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else {
    #trim vector elements to remove leading and trailing blanks
    categ.init <- trim(categ.init)
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of 
  # raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for (k in 1:q) {
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      } else
        agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # creating the rxq rater-category matrix representing the distribution of 
  # subjects by rater and category
  classif.mat <- matrix(0,nrow=r,ncol=q)
  for (k in 1:q) {
    if (is.numeric(ratings.mat)) {
      with.mis <-(t(ratings.mat)==categ[k])
      without.mis <- replace(with.mis,is.na(with.mis),FALSE)
      classif.mat[,k] <- without.mis%*%rep(1,n)
    } else
      classif.mat[,k] <- (t(ratings.mat)==categ[k])%*%rep(1,n)
  }
  
  # calculating conger's kappa coefficient
  ri.vec <- agree.mat %*% rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1)) %*% rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  ng.vec <- classif.mat %*% rep(1,q)
  pgk.mat <- classif.mat/(ng.vec %*% rep(1,q))
  p.mean.k <- (t(pgk.mat) %*% rep(1,r))/r 
  s2kl.mat <- (t(pgk.mat) %*% pgk.mat - r * p.mean.k %*% t(p.mean.k))/(r-1)
  pe <- sum(weights.mat * (p.mean.k%*%t(p.mean.k) -  s2kl.mat/r))
  conger.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of conger's kappa coefficient
  bkl.mat <- (weights.mat+t(weights.mat))/2
  pe.ivec1 <- r*(agree.mat%*%t(t(p.mean.k)%*%bkl.mat))
  pe.ivec2 = rep(0,n)
  if (is.numeric(ratings.mat)) {
    for(l in 1:q) {
      l.mis <-(ratings.mat==categ[l])
      delta.ig.mat <- replace(l.mis,is.na(l.mis),FALSE)
      pe.ivec2 <- pe.ivec2 + delta.ig.mat%*%(pgk.mat%*%bkl.mat[,l]) 
    }
  } else {
    for (k in 1:q) {
      delta.ig.mat <- (ratings.mat==categ[k])
      pe.ivec2 <- pe.ivec2 + delta.ig.mat%*%(pgk.mat%*%bkl.mat[,l])
    }
  }
  
  pe.ivec <- (pe.ivec1-pe.ivec2)/(r*(r-1)) 
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next 
  # ratio calculation always possible.
  den.ivec <- den.ivec - (den.ivec==0)
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  conger.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe) 
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.conger <- ((1-f)/(n*(n-1))) * sum((conger.ivec.x - conger.kappa)^2)
  # conger's kappa standard error
  stderr <- sqrt(var.conger)
  p.value <- 2*(1-pt(conger.kappa/stderr,n-1))
  
  # confidence bounds
  lcb <- conger.kappa - stderr*qt(1-(1-conflev)/2,n-1)
  ucb <- min(1,conger.kappa + stderr*qt(1-(1-conflev)/2,n-1))
  
  if (print) {
    cat("Conger's Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement: ',pa,'Percent chance agreement: ',pe,'\n')
    cat("Conger's kappa coefficient: ",conger.kappa,'Standard error:',
        stderr,'\n')
    if (weights!="unweighted") {
      if (!is.numeric(weights)) {
        cat('Weights: ', weights,'\n')
      } else
        cat('Weights: Custom Weights\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(c(pa,pe,conger.kappa,stderr,p.value))
}

#==============================================================================
#bp.coeff.raw: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and 
# its standard error for multiple raters when input dataset is a nxr matrix 
# of alphanumeric ratings from n subjects and r raters 
#-------------
# The input data "ratings" is a nxr data frame of raw alphanumeric ratings
# from n subjects and r raters. Exclude all subjects that are not rated by any 
# rater.
# 
# Bibliography:
# Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, 
#     misuses, and alternatives." Educational and Psychological Measurement, 41,
#     687-699.
#==============================================================================

#' Brennan-Prediger kappa coefficient
#' 
#' Computes Brennan-Prediger kappa coefficient and standard error for multiple
#' raters when data is an \code{n x r} matrix of alphanumeric ratings from 
#' \code{n} subjects and \code{r} raters, excluding all subjects that are not 
#' rated by any rater.
#' 
#' \code{weight} is an option matrix of weights or one of \code{"quadratic"},
#' \code{"linear"}, \code{"ordinal"}, \code{"radical"}, \code{"ratio"},
#' \code{"circular"}, \code{"bipolar"}, or \code{"unweighted"} for an identity
#' matrix.
#' 
#' @param ratings an \code{n x r} matrix of raw alphanumeric ratings from
#' \code{n} subjects and \code{r} raters
#' @param weights optional weighting; see details
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

bp.coeff.raw <- function(ratings, weights = "unweighted",
                         conflev = 0.95, N = Inf, print = TRUE) {
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N               # final population correction 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else {
    # trim vector elements to remove leading and trailing blanks
    categ.init <- trim(categ.init)
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of 
  # raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)) {
    k.mis <-(ratings.mat==categ[k])
    in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
	} else
    agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient
  ri.vec <- agree.mat %*% rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1)) %*% rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n)) %*% (agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  den.ivec <- ri.vec*(ri.vec-1)
  # this operation replaces each 0 value with -1 to make the next ratio 
  # calculation always possible.
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
    
    if (weights!="unweighted") {
      if (!is.numeric(weights)) {
        cat('Weights: ', weights,'\n')
      } else
        cat('Weights: Custom Weights\n')
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
