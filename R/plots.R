### plots
# ba_plot
###

#' Bland-Altman plot
#' 
#' Creates a scatter plot using the method of Bland and Altman where the x-
#' and y-axes are the means of and differences between \code{x} and \code{y},
#' respectively. Lines for the mean and \code{+/-1.96} times the standard
#' deviation of the difference between \code{x} and \code{y} are shown.
#' 
#' @param x,y measurements
#' @param trend logical; if \code{TRUE}, adds trend line to data
#' @param span numeric value that controls the degree of smoothing; see
#' \code{\link{loess}}
#' @param trend.n integer value for the number of data points to be used by
#' \code{\link{predict.loess}}
#' @param xlab,ylab x- and y-axis labels
#' @param xlim,ylim x- and y-axis limits
#' @param panel.first an expression to be evaluated after the plot axes are
#' set up but before any plotting takes place (default is just a \code{grid});
#' see \code{\link{plot.default}}
#' @param ... additional parameters passed to \code{\link{plot}}
#' @param line_pars additional graphical parameters passed to \code{link{text}}
#' or \code{\link{abline}}; all parameters will be passed as-is except
#' \code{col} which should be a vector of length one, two, or three to color
#' the lines and text accordingly; see examples
#' 
#' @references
#' Altman DG, Bland JM (1983). "Measurement in medicine: the analysis of method
#' comparison studies". \emph{The Statistician} \strong{32}: 307–317.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100, .5)
#' 
#' ## basic usage
#' ba_plot(x, y)
#' 
#' ## use line_pars to set graphical parameters of individual lines
#' ba_plot(x, y, col = c('black','green')[(abs(x - y) > 2) + 1],
#'         panel.first = NULL, pch = '.', cex = 5,
#'         line_pars = list(col = c('yellow', 'red')))
#'         
#' ## panel.first/panel.last is used as in plot.default
#' ba_plot(x, y, cex = 2,
#'         col = as.numeric(cut(x - y, c(-Inf, -2, 2, Inf))) + 2,
#'         panel.first = {
#'           p <- par('usr')
#'           rect(p[1], p[3], p[2], p[4], col = 'cornsilk')
#'           grid(col = 'grey90', lty = 1)
#'         },
#'         line_pars = list(col = 3:5, family = 'HersheySerif'))
#' 
#' @export

ba_plot <- function(x, y, trend = TRUE, span = 0.5, trend.n = 1000L,
                    xlab = '', ylab = '', xlim = NULL, ylim = NULL,
                    panel.first = grid(), panel.last = NULL, ..., 
                    line_pars = list(col = 1:2, lty = c(2,1,2))) {
  xl <- as.character(substitute(x))
  yl <- as.character(substitute(y))
  dd <- data.frame(x = (x + y) / 2, y = y - x)
  cols <- line_pars$col
  cols <- if (is.null(cols))
    c('dodgerblue2','firebrick','dodgerblue2') else rep_len(cols, 3)
  line_pars$col <- NULL
  
  ## text positions
  nums <- sort(c(sd(dd$y) * 1.96 * c(-1, 1) + mean(dd$y), mean(dd$y)))
  txt <- diff(range(dd$x)) * .1
  
  plot(dd, ..., xlab = xlab, ylab = ylab,
       xlim = xlim %||% range(dd$x) + c(diff(nums)) * c(0, .15),
       ylim = ylim %||% range(dd$y) + c(diff(nums)) * c(-.1, .1),
       panel.first = {
         panel.first
         p <- par('usr')
         do.call('abline', c(list(h = nums, col = cols), line_pars))
         do.call('text', c(list(
           x = p[2] - txt, y = nums, col = cols, pos = 3,
           labels = c('-1.96 SD', 'mean', '+1.96 SD')), line_pars))
         do.call('text', c(list(
           x = p[2] - txt, y = nums, col = cols, pos = 1,
           labels = roundr(nums, 2)), line_pars))
       },
       panel.last = {
         panel.last
         if (trend) {
           lo <- loess(y ~ x, dd, span = span, degree = 1)
           sx <- seq(min(dd$x), max(dd$x), length.out = trend.n)
           lines(sx, predict(lo, sx), col = 'red', lwd = 2)
         }
       })
  title(xlab = xlab %||% sprintf('Average of %s and %s', xl, yl),
        ylab = ylab %||% sprintf('Difference between %s and %s', xl, yl))
  invisible(dd)
}
