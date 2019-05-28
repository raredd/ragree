### agree stats and utilities
# gkappa, unalike, unalike.default, unalike.matrix, unalike.data.frame,
# maxwell, maxwelln
# 
# unexported:
# unalike1, unalike2, untable
###


#' Generalized kappa
#' 
#' Calculates generalized kappa for k-raters
#' 
#' @param ratings \code{n * m} matrix or data frame; \code{n} subjects, 
#' \code{m} raters
#' @param data data in long format with three columns: ids, raters, scores
#' @param id,rater,score column names corresponding to IDs, raters, and scores
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
#' \tab \code{$detail} \tab  \cr
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
#' library('irr')
#' data(diagnoses)
#' 
#' gkappa(ratings = diagnoses)
#' 
#' ## fleiss's variance calculation is slightly different
#' kappam.fleiss(diagnoses)
#' 
#' @export

gkappa <- function(ratings, data, id = 'id', rater = 'rater', score = 'score') {
  
  if (!missing(data)) {
    K   <- length(unique(data[[rater]]))  ## raters
    R   <- length(unique(data[[score]]))  ## categories
    N   <- nrow(data) / K                 ## patients
    tbl <- ftable(data[[id]], data[[score]])[seq.int(N), , drop = FALSE]
  }
  
  if (!missing(ratings)) {
    ratings <- as.matrix(ratings)
    N   <- nrow(ratings)
    K   <- ncol(ratings)
    R   <- length(unique(c(ratings)))
    tbl <- ftable(data.frame(id    = rep(seq.int(N), K),
                             score = c(ratings)))[seq.int(N), , drop = FALSE]
  }
  
  ## proportion of all classifications assigned to category j
  p_hat_j <- function(...) (colSums(tbl) / (N * K))[c(...)]
  
  ## proportion of pairs agreeing for each person i 
  p_hat_i <- rowSums(tbl  * (tbl - 1) / (K * (K - 1)))
  
  ## overall agreement
  p_bar <- sum(p_hat_i) / N
  
  ## expected chance agreement
  p_bar_e <- sum(p_hat_j(1:R) ** 2)
  
  ## generalized kappa statistic
  K_hat_G <- (p_bar - p_bar_e) / (1 - p_bar_e)
  
  ## variance estimate
  sigma2_Kg <- (2 / (N * K * (K - 1))) *
    (sum(p_hat_j(1:R) ** 2) - (2 * K - 3) *
       (sum(p_hat_j(1:R) ** 2) ** 2) + 2 *
       (K - 2) * sum(p_hat_j(1:R) ** 3)) /
    (1 - sum(p_hat_j(1:R) ** 2)) ** 2
  
  ## test statistic
  Z <- K_hat_G / sqrt(sigma2_Kg)
  
  ## comparing test statistic to a standard normal distribution
  pval <- 2 * (1 - pnorm(abs(Z)))
  
  res <- list(
    method      = sprintf('Generalized kappa for %s raters', K),
    ragree.name = 'Kappa',
    subjects    = N,
    raters      = K,
    categories  = R,
    value       = K_hat_G,
    variance    = sigma2_Kg,
    statistic   = Z,
    stat.name   = 'z',
    p.value     = pval, 
    detail      = list(p_hat_j = p_hat_j(seq.int(R)),
                       p_bar = p_bar,
                       p_bar_e = p_bar_e)
  )
  
  structure(res, class = c('ragree', 'list'))
}

#' Coefficient of unalikeability
#' 
#' Function to calculate the unalikeability coefficient to quantify the amount 
#' of variability in categorical data.
#' 
#' @param x a vector of categorical data
#' 
#' alternatively, an \code{n x r} matrix of \code{n} subjects by \code{r}
#' raters or an \code{n x m} data frame of \code{n} subjects; when \code{x}
#' is a data  frame, the user must specify which columns correspond to
#' \code{id}, \code{rater}, and \code{score}
#' @param ... additional arguments passed to or from other methods
#' @param method the method for calculating the unalikeability coefficient;
#' see details
#' @param id,rater,score column names corresponding to IDs, raters, and scores
#' @param summary logical; if \code{TRUE}, prints summary statistics for the
#' unalikeability coefficients
#' @param plot logical; if \code{TRUE}, prints a heat map of unalikeability
#' coefficients
#' 
#' @return A list containing the following:
#' 
#' \tabular{rlll}{
#' \tab \code{$method} \tab agreement method \cr
#' \tab \code{$ragree.name} \tab method type \cr
#' \tab \code{$subjects} \tab number of subjects \cr
#' \tab \code{$raters} \tab number of raters \cr
#' \tab \code{$categories} \tab number of categories \cr
#' \tab \code{$value} \tab median of all unalikeability coefficients \cr
#' \tab \code{$summary} \tab a data frame with the summary information
#' printed when \code{summary = TRUE} \cr
#' \tab \code{$data} \tab a long data frame with all coefficients \cr
#' }
#' 
#' @details
#' The coefficient of unalikeability describes a concept of variability for 
#' categorical variables and provides a quantitative method for its 
#' measurement. A smaller coefficient is better corresponding to less variation
#' in the scores.
#' 
#' For the case of a finite number of observations (\code{n}), a finite number 
#' of categories (\code{m}) and a finite number of objects, \eqn{k_i}, within
#' category \code{i}, will allow expression of the coefficient of unalikeablity
#' as:
#' 
#' \deqn{u = 1 - \sum p_i ^ 2} where \eqn{p_i = k_i / n}.
#' 
#' The interpretation of \eqn{u} is that it represents the proportion of 
#' possible comparisons (pairings) which are unalike. Note that \eqn{u}
#' includes comparisons of each response with itself.
#' 
#' Currently, two methods for calculating the coefficient are implemented. If
#' \code{method = 1}, then the formula described above is used. If
#' \code{method = 2}, then the formula described in Perry (2005).
#' 
#' @seealso
#' \code{RcmdrPlugin.ISCSS::unalike}
#' 
#' @author
#' Robert Redd \email{rredd@@jimmy.harvard.edu}
#' 
#' @references
#' Kader, GD. Variability for Categorical Variables. \emph{Journal of
#' Statistics Education}, Vol. 15, No. 2 (2007).
#' 
#' Perry, M. and Kader, G. Variation as Unalikeability. \emph{Teaching
#' Statistics}, Vol. 27, No. 2 (2005), pp. 58-60.
#' 
#' @examples
#' unalike(1, 2)
#' unalike(rep(1, 10))
#' 
#' 
#' ## examples in Kader (2007):
#' l <- list(
#'   group1 = rep(c('A', 'B'), c(7, 3)),
#'   group2 = rep(c('A', 'B'), c(5, 5)),
#'   group3 = rep(c('A', 'B'), c(1, 9)),
#'   group4 = rep(c('A', 'B', 'C'), c(2, 3, 5))
#' )
#' 
#' sapply(l, unalike)
#' 
#' 
#' ## matrix/data frames are assumed to be subjects x raters
#' mat <- do.call('cbind', l[1:3])
#' unalike(mat) ## see Kader
#' 
#' 
#' library('irr')
#' data(diagnoses)
#' 
#' kappam.fleiss(diagnoses)
#' unalike(as.matrix(diagnoses))
#' 
#' library('ggplot2')
#' unalike(as.matrix(diagnoses), plot = TRUE)
#' 
#' 
#' dat <- data.frame(
#'   id    = rep(seq.int(nrow(diagnoses)), ncol(diagnoses)),
#'   rater = rep(names(diagnoses), each = nrow(diagnoses)),
#'   score = unlist(diagnoses)
#' )
#' unalike(dat)
#' 
#' @export

unalike <- function(x, ...) {
  UseMethod('unalike')
}

unalike1 <- function(x) {
  ## helper function to calculate the unalikeability coefficient
  x <- if (inherits(x, 'table'))
    x else table(x)
  1 - sum(prop.table(x) ^ 2)
}

unalike2 <- function(x) {
  x <- if (inherits(x, 'table'))
    untable(x) else x
  n <- length(x)
  o <- outer(x, x, `!=`)
  # mean(o[lower.tri(o)])
  sum(o[row(o) != col(o)]) / (n ^ 2 - n)
}

untable <- function(x) {
  stopifnot(
    inherits(x, 'table'),
    length(dim(x)) == 1L
  )
  rep(seq_along(x), x)
}

#' @rdname unalike
#' @export
unalike.default <- function(x, ..., method = 1L) {
  x <- if (inherits(x, 'table'))
    x else c(x, ...)
  
  if (method == 1L)
    unalike1(x)
  else if (method == 2L)
    unalike2(x)
  else stop('Invalid method - should be 1 or 2', call. = FALSE)
}

#' @rdname unalike
#' @export
unalike.matrix <- function(x, ...) {
  x <- as.data.frame(cbind(`_id_` = seq.int(nrow(x)), x))
  x <- reshape(x, idvar = '_id_', varying = list(2:ncol(x)),
               direction = 'long')
  names(x) <- c('id', 'rater', 'score')
  
  unalike(x, id = 'id', rater = 'rater', score = 'score', ...)
}

#' @rdname unalike
#' @export
unalike.data.frame <- function(x, id = 'id', rater = 'rater', score = 'score',
                               summary = TRUE, plot = FALSE, ...) {
  if (any(idx <- !(c(id, rater, score) %in% names(x))))
    stop(
      sprintf('%s not found in data', toString(c(id, rater, score)[idx]))
    )
  
  ## calculate the unalikeability coefficient for each id
  x$unalike <- ave(as.numeric(as.factor(x[[score]])), list(x[[id]]),
                   FUN = function(ii) unalike(ii, ...))
  ua <- x[!duplicated(x[[id]]), ]$unalike
  
  res <- data.frame(
    LCI    = quantile(ua, probs = 0.025),
    Min    = min(ua),
    Median = median(ua),
    Mean   = mean(ua),
    Max    = max(ua),
    UCI    = quantile(ua, probs = 0.975)
  )
  
  if (summary) {
    cat('\nSummary of unalikeability coefficients:\n\n')
    print(res, row.names = FALSE, digits = 3L)
    cat('\n\n')
  }
  
  if (plot) {
    x$id <- as.factor(x$id)
    
    p <- ggplot(x, aes(x = factor(rater), y = rev(id))) + 
      geom_tile(aes(fill = unalike), colour = 'white') + 
      scale_fill_gradient(name = 'Un-alikeability\ncoefficient',
                          limits = c(0, 1), 
                          low = 'white', 
                          high = 'steelblue') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(labels = rev(unique(x$id)), expand = c(0, 0)) +
      theme_bw() + 
      theme(legend.position = 'right',
            axis.text.x = element_text(size = 15, colour = 'grey50'),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    print(p)
  }
  
  res <- list(
    method      = sprintf('Unalikeability coefficient for %s raters',
                          nrow(x) / length(unique(x[[id]]))),
    ragree.name = 'Unalike',
    subjects    = length(unique(x[[id]])),
    raters      = nrow(x) / length(unique(x[[id]])),
    categories  = length(unique(x[[score]])),
    value       = median(x[['unalike']]),
    summary     = res,
    data        = x
  )
  
  structure(res, class = 'ragree')
}

#' Maxwell's RE
#' 
#' Computes Maxwell's RE coefficient for binary data of two raters.
#' 
#' @param ratings an \code{n x 2} matrix or data frame with \code{n} subjects
#' and 2 raters with binary ratings
#' 
#' note if \code{ratings} contains more than two columns, the first two are
#' used and the remaining are ignored; to compute the pairwise coefficient for
#' an \code{n x r} matrix with \code{r} raters, use \code{maxwelln}; see
#' examples
#' 
#' @references
#' Maxwell, A.E. (1977). Coefficients of agreement between observers and their
#' interpretation. \emph{Br J Psychiatry} \strong{130}: 79-83.
#' 
#' @seealso
#' \code{\link[irr]{maxwell}}
#' 
#' @examples
#' anx <- +(anxiety > 1)
#' maxwell(anx)
#' 
#' ## compare
#' irr::maxwell(anx[, -3])
#' 
#' ## to get RE for all pairs of raters
#' maxwelln(anx)
#' 
#' @export

maxwell <- function(ratings) {
  x <- na.omit(as.matrix(ratings))
  n <- nrow(x)
  
  r1 <- as.factor(x[, 1L])
  r2 <- as.factor(x[, 2L])
  
  lvl <- unique(c(levels(r1), levels(r2)))
  tbl <- table(factor(r1, lvl), factor(r2, lvl))
  
  if (!identical(dim(tbl), c(2L, 2L)))
    stop('ratings must be binary', call. = FALSE)
  
  list(
    method = 'Maxwell\'s RE for 2 raters',
    ragree.name = 'RE',
    subjects = n,
    raters = 2L,
    categories = 2L,
    value = 2 * sum(diag(tbl)) / n - 1
  )
}

#' @rdname maxwell
#' @export
maxwelln <- function(ratings) {
  nr <- ncol(ratings)
  cc <- combn(nr, 2L)
  
  res <- lapply(seq.int(ncol(cc)), function(ii) {
    maxwell(ratings[, cc[, ii]])
  })
  
  names(res) <- apply(cc, 2L, function(ii)
    paste0(colnames(ratings)[ii], collapse = '.vs.'))
  
  do.call('cbind', res)
}
