#' Generalized kappa
#' 
#' Calculates generalized kappa for k-raters
#' 
#' @usage 
#' gkappa(ratings, data, id = 'id', rater = 'rater', score = 'score')
#' 
#' @param ratings \code{n * m} matrix or data frame; \code{n} subjects, 
#' \code{m} raters
#' @param data data in long format with three columns: ids, raters, scores
#' @param id name of the column with ids
#' @param rater name of the column with raters
#' @param score name of the column with scores
#' 
#' @return A list containing the generalized kappa coefficient, the estimated 
#' variance, the test statistic, and corresponding p-value comparing the 
#' statistic to a standard normal distribution.
#' 
#' \tabular{rlll}{
#' \tab \code{$method} \tab agreement method \cr
#' \tab \code{$ragree.name} \tab method type \cr
#' \tab \code{$subjects} \tab number of subjects \cr
#' \tab \code{$raters} \tab number of raters \cr
#' \tab \code{$categories} \tab number of categories \cr
#' \tab \code{$value} \tab generalized kappa statistic \cr
#' \tab \code{$variance} \tab variance \cr
#' \tab \code{$statistic} \tab test statistic \cr
#' \tab \code{$stat.name} \tab test statistic name \cr
#' \tab \code{$p.value} \tab p-value comparing statistic to standard normal 
#' distribution \cr \cr
#' \tab \code{\strong{$detail}} \tab \strong{a list of three containing:} \cr
#' \tab \code{...$p_hat_j} \tab proportion of all classifications assigned to 
#' category j; the pairwise agreement \cr
#' \tab \code{...$p_bar} \tab overall agreement \cr
#' \tab \code{...$p_bar_e} \tab expected chance agreement \cr
#' }
#' @author Robert Redd \email{rredd@@jimmy.harvard.edu}
#' 
#' @examples
#' ## long format using data argument
#' gkappa(data = disease.l, id = 'patient', rater = 'rater', score = 'category')
#' 
#' ## wide format using ratings argument
#' gkappa(ratings = disease.w)
#' 
#' ## example from irr package
#' library(irr)
#' data(diagnoses)
#' 
#' gkappa(ratings = diagnoses)
#' 
#' ## fleiss's variance calculation is slightly different
#' kappam.fleiss(diagnoses)
#' @export

gkappa <- function(ratings, data, id = 'id', rater = 'rater', score = 'score') {
  
  if (!missing(data)) {
    K <- length(unique(data[[rater]]))  ## raters
    R <- length(unique(data[[score]]))  ## categories
    N <- nrow(data) / K                 ## patients
    tmp <- ftable(data[[id]], data[[score]])[1:N, , drop = FALSE]
  }
  
  if (!missing(ratings)) {
    N <- nrow(ratings)
    K <- ncol(ratings)
    R <- length(unique(unlist(ratings)))
    tmp <- ftable(data.frame(id = rep(1:N, K),
                             score = unlist(ratings)))[1:N, , drop = FALSE]
  }
  
  ## proportion of all classifications assigned to category j
  p_hat_j <- function(...) (colSums(tmp) / (N * K))[c(...)]
  
  ## proportion of pairs agreeing for each person i 
  p_hat_i <- rowSums(tmp  * (tmp - 1) / (K * (K - 1)))
  
  ## overall agreement
  p_bar <- sum(p_hat_i) / N

  ## expected chance agreement
  p_bar_e <- sum(p_hat_j(1:R) ** 2)
  
  ## generalized kappa statistic
  K_hat_G <- (p_bar - p_bar_e) / (1 - p_bar_e)
  
  ## variance estimate
  sigma2_Kg <-   (2 / (N * K * (K - 1))) * 
    (sum(p_hat_j(1:R) ** 2) - (2 * K - 3) * 
       (sum(p_hat_j(1:R) ** 2) ** 2) + 2 * 
       (K - 2) * sum(p_hat_j(1:R) ** 3)) / 
    (1 - sum(p_hat_j(1:R) ** 2)) ** 2
  
  ## test statistic
  Z <- K_hat_G / sqrt(sigma2_Kg)
  
  ## comparing test statistic to a standard normal distribution
  pval <- 2 * (1 - pnorm(abs(Z)))
  
  zzz <- list(method = 'Generalized kappa for m raters',
              ragree.name = 'Kappa',
              subjects = N,
              raters = K,
              categories = R,
              value = K_hat_G,
              variance = sigma2_Kg,
              statistic = Z,
              stat.name = 'z',
              p.value = pval, 
              detail = list(p_hat_j = p_hat_j(1:R),
                            p_bar = p_bar,
                            p_bar_e = p_bar_e))
  class(zzz) <- c('ragree', 'list')
  
  return(zzz)
}
