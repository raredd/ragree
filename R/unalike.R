#' Coefficient of unalikeability
#' 
#' Function to calculate the unalikeability coefficient to quantify the amount 
#' of variability in categorical data.
#' 
#' @usage
#' unalike(...)
#' 
#' ## S3 method for class 'matrix':
#' unalike(x, ...)
#' 
#' ## S3 method for class 'data.frame':
#' unalike(x, id, rater, score, summary = TRUE, plot = FALSE)
#' 
#' @param ... numeric values or vector representing the count of each variable
#' @param x an \code{n x r} matrix of \code{n} subjects by \code{r} raters or
#' an \code{n x m} data frame of \code{n} subjects; when \code{x} is a data 
#' frame, the user must specify which columns correspond to \code{id}, 
#' \code{rater}, and \code{score}
#' @param id data frame column name corresponding to IDs
#' @param rater data frame column name corresponding to raters
#' @param score data frame column name corresponding to rater scores
#' @param summary logical; if \code{TRUE}, prints summary statistics for the
#' unalikeability coefficients
#' @param plot logical; if \code{TRUE}, prints a heat map of unalikeability
#' coefficients; see details
#' 
#' @return A list containing the following:
#' 
#' \tabular{rlll}{
#' \tab \code{$method} \tab agreement method \cr
#' \tab \code{$ragree.name} \tab method type \cr
#' \tab \code{$subjects} \tab number of subjects \cr
#' \tab \code{$raters} \tab number of raters \cr
#' \tab \code{$categories} \tab number of categories \cr
#' \tab \code{$value} \tab mean of all unalikeability coefficients \cr
#' \tab \code{$zzz} \tab a data frame with the summary information printed
#' when \code{summary = TRUE} \cr
#' \tab \code{$dat} \tab a long data frame with all coefficients \cr
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
#' If \code{plot} is \code{TRUE}, a \code{\link[ggplot2]{ggplot}} heatmap will
#' be printed.
#'
#' @author Robert Redd \email{rredd@@jimmy.harvard.edu}
#' @references Kader, GD. Variability for Categorical Variables. \emph{Journal 
#' of Statistics Education}, Vol. 15, No. 2: 2007.
#' 
#' @examples
#' ## examples in Kader:
#' ## sample data
#' l <- list(grp1 = c(rep('A', 7), rep('B', 3)),
#'           grp2 = rep(c('A','B'), each = 5),
#'           grp3 = c(rep('A', 1), rep('B', 9)))
#' 
#' unalike(1, 9)
#' sapply(l, function(x) unalike(x))
#' unalike(sample(c('a','b','c'), 20, replace = TRUE))
#' 
#' mat <- do.call(cbind, l)
#' colnames(mat) <- c('g1','g2','g3')
#' unalike(t(mat)) ## see Kader
#' 
#' \dontrun{
#' library(irr)
#' data(diagnoses)
#' kappam.fleiss(diagnoses)
#' unalike(as.matrix(diagnoses), plot = TRUE)
#' }
#' @export

unalike <- function(...) UseMethod('unalike')

#' @export
unalike.default <- function(...) {
  input <- list(...)
  
  ## helper function to calculate the unalikeability coefficient
  unalike.helper <- function(...)
    round(1 - sum(((c(...)) / sum(c(...))) ** 2), 3)
  
  ## if there are multiple args in ...
  ## if a list of numbers is passed to ...
  if (all(class(unlist(input)) %in% c('numeric', 'integer')))
    return(unalike.helper(...))
  ## if a list of categorical variables are passed
  if (is.character(unlist(input))) {
    warning('character values coerced to numeric')
    return(unalike.helper(table(as.numeric(as.factor(...)))))
  }
  ## if a list of factors are passed
  if (is.factor(unlist(input))) {
    warning('factor levels coerced to numeric')
    return(unalike.helper(table(as.numeric(...))))
  }
}

#' @export
unalike.matrix <- function(x, ...) {
  x <- as.data.frame(cbind(id = 1:nrow(x), x))
  x <- reshape(x, idvar = 'id', varying = list(2:ncol(x)), direction = 'long')
  names(x) <- c('id','rater','score')
  unalike.data.frame(x, id = 'id', rater = 'rater', score = 'score', ...)
}

#' @export
unalike.data.frame <- function(x, id, rater, score, 
                               summary = TRUE, plot = FALSE) {

  if (any(idx <- !(c(id, rater, score) %in% names(x))))
    stop(sprintf('%s not found in data',
                 paste(c(id, rater, score)[idx], collapse = ', ')))
  
  ## helper function to calculate the unalikeability coefficient
  unalike.helper <- function(...)
    round(1 - sum(((c(...)) / sum(c(...))) ** 2), 3)
  
  ## calculate the unalikeability coefficient for each id
  x$unalike <- ave(as.numeric(as.factor(x[[score]])), list(x[[id]]),
                   FUN = function(ii) unalike.helper(table(ii)))
  zzz1 <- x[!duplicated(x[[id]]), ]$unalike
  
  zzz <- data.frame(Minimum = min(zzz1),
                    Mean = mean(zzz1),
                    Median = median(zzz1),
                    Maximum = max(zzz1),
                    '95% CI' = sprintf('(%s, %s)', 
                                       round(quantile(zzz1, probs = .025), 3),
                                       round(quantile(zzz1, probs = .975), 3)),
                    check.names = FALSE, stringsAsFactors = FALSE)
  
  if (summary) {
    cat('\nSummary of unalikeability coefficients:\n\n')
    print(zzz, row.names = FALSE)
    cat('\n\n')
  }
  
  if (plot) {
    require(ggplot2)
    
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
  
  zzz <- list(method = 'Unalikeability coefficient for m raters',
              ragree.name = 'Unalike',
              subjects = length(unique(x[[id]])),
              raters = nrow(x) / length(unique(x[[id]])),
              categories = length(unique(x[[score]])),
              value = mean(x[['unalike']]),
              zzz = zzz,
              dat = x)
  class(zzz) <- c('ragree', 'list')
  
  return(zzz)
}
