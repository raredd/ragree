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
#' Bland, Altman. Sometime.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100, .5)
#' ba_plot(x, y, pch = 19)
#' 
#' ba_plot(x, y, col = c('black','green')[(abs(x - y) > 2) + 1],
#'         panel.first = NULL, pch = '.', cex = 5,
#'         line_pars = list(col = c('yellow', 'red')))
#'         
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

ba_plot <- function(x, y, xlab = '', ylab = '', xlim = NULL, ylim = NULL,
                    panel.first = grid(), ..., 
                    line_pars = list(col = c('dodgerblue2', 'firebrick'),
                                     lty = 'dashed', lwd = 2, cex = 1.5)) {
  xl <- as.character(substitute(x))
  yl <- as.character(substitute(y))
  dd <- data.frame(x = (x + y) / 2, y = x - y)
  cols <- line_pars$col
  cols <- if (is.null(cols))
    c('dodgerblue2','firebrick','dodgerblue2') else recycle(1:3, cols)
  line_pars$col <- NULL
  
  ## text positions
  nums <- sort(c(sd(dd$y) * 1.96 * c(-1, 1) + mean(dd$y), mean(dd$y)))
  txt <- diff(range(dd$x)) * .1
  
  plot(dd, ..., xlab = xlab, ylab = ylab,
       xlim = xlim %||% range(dd$x) + c(diff(nums)) * c(0, .15),
       ylim = ylim %||% range(dd$y) + c(diff(nums)) * c(-.1, .1),
       panel.first = {
         p <- par('usr')
         panel.first
         do.call('abline', c(list(h = nums, col = cols), line_pars))
         do.call('text', c(list(
           x = p[2] - txt, y = nums, col = cols, pos = 3,
           labels = c('-1.96 SD','mean','+1.96 SD')), line_pars))
         do.call('text', c(list(
           x = p[2] - txt, y = nums, col = cols, pos = 1,
           labels = roundr(nums, 2)), line_pars))
       })
  title(xlab = xlab %||% sprintf('Average of %s and %s', xl, yl),
        ylab = ylab %||% sprintf('Difference between %s and %s', xl, yl))
  invisible(dd)
}
