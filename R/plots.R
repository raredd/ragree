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
#' @param trend character string giving the type of trend line, one of
#' \code{"loess"} (default), \code{"gam"} (requires the \pkg{\link{mgcv}}
#' package), or \code{"none"} where no curve is drawn
#' @param span if \code{trend = "loess"}, a numeric value that controls the
#' degree of smoothing (see \code{\link{loess}})
#' @param trend.n integer value for the number of data points to be used by
#' \code{\link{predict.loess}}
#' @param alpha confidence level affecting the critical values and trend
#' confidence bands
#' @param xlab,ylab x- and y-axis labels
#' @param xlim,ylim x- and y-axis limits
#' @param panel.first an expression to be evaluated after the plot axes are
#' set up but before any plotting takes place (default is just a \code{grid});
#' see \code{\link{plot.default}}
#' @param panel.last an expression to be evaluated after plotting has taken
#' place but before the axes, title, and box are added; see comments about
#' \code{panel.first}
#' @param line_pars additional graphical parameters passed to \code{\link{text}}
#' or \code{\link{abline}}; all parameters will be passed as-is except
#' \code{col} which should be a vector of length one, two, or three to color
#' the lines and text accordingly; see examples
#' @param ... additional parameters passed to \code{\link{plot}}
#' 
#' @references
#' Altman DG, Bland JM (1983). Measurement in medicine: the analysis of method
#' comparison studies. \emph{The Statistician} \strong{32}: 307-317.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100, .5)
#' 
#' ## basic usage
#' ba_plot(x, y)
#' ba_plot(x, y, trend = 'gam')
#' ba_plot(x, y, trend = 'none', alpha = 0.1)
#' 
#' 
#' ## use line_pars to set graphical parameters of individual lines
#' ba_plot(x, y, col = c('black', 'green')[(abs(x - y) > 2) + 1L],
#'         panel.first = NULL, pch = '.', cex = 5,
#'         line_pars = list(col = c('yellow', 'red')))
#'         
#' ## panel.first/panel.last is used as in plot.default
#' ba_plot(x, y, cex = 2,
#'         col = as.numeric(cut(x - y, c(-Inf, -2, 2, Inf))) + 2,
#'         panel.first = {
#'           p <- par('usr')
#'           rect(p[1L], p[3L], p[2L], p[4L], col = 'cornsilk')
#'           grid(col = 'grey90', lty = 1)
#'         },
#'         line_pars = list(col = 3:5, family = 'HersheySerif'))
#' 
#' @export

ba_plot <- function(x, y, trend = c('loess', 'gam', 'none'),
                    span = 0.5, trend.n = 1000L, alpha = 0.05,
                    xlab = '', ylab = '', xlim = NULL, ylim = NULL,
                    panel.first = grid(), panel.last = NULL, ..., 
                    line_pars = list(col = 1:2, lty = c(2,1,2))) {
  trend <- match.arg(trend)
  cv <- qnorm(alpha / 2, lower.tail = FALSE)
  cc <- paste0(c('-', '+'), sprintf('%.2f SD', cv))
  
  xl <- as.character(substitute(x))
  yl <- as.character(substitute(y))
  dd <- data.frame(x = (x + y) / 2, y = y - x)
  
  cols <- line_pars$col
  cols <- if (is.null(cols))
    c('dodgerblue2', 'firebrick', 'dodgerblue2') else rep_len(cols, 3L)
  line_pars$col <- NULL
  
  ## text positions
  nums <- sort(c(sd(dd$y) * cv * c(-1, 1) + mean(dd$y), mean(dd$y)))
  txt  <- diff(range(dd$x)) * .1
  
  plot(dd, ..., xlab = xlab, ylab = ylab,
       xlim = xlim %||% range(dd$x) + c(diff(nums)) * c(0, .15),
       ylim = ylim %||% range(dd$y) + c(diff(nums)) * c(-.1, .1),
       
       panel.first = {
         panel.first
         p <- par('usr')
         do.call('abline', c(list(h = nums, col = cols), line_pars))
         do.call('text', c(list(
           x = p[2L] - txt, y = nums, col = cols, pos = 3,
           labels = c(cc[1L], 'mean', cc[2L])), line_pars))
         do.call('text', c(list(
           x = p[2L] - txt, y = nums, col = cols, pos = 1,
           labels = roundr(nums, 2L)), line_pars))
       },
       
       panel.last = {
         if (trend != 'none') {
           ## smoothing
           nx <- seq(min(dd$x), max(dd$x), length.out = trend.n)
           
           do_poly <- function(x, y1, y2) {
             polygon(c(x, rev(x)), c(y1, rev(y2)), border = NA,
                     col = adjustcolor(2, alpha.f = 0.2))
           }
           
           if (trend == 'loess') {
             lo <- loess(y ~ x, dd, span = span, degree = 1)
             pr <- predict(lo, nx, se = TRUE)
             lines(nx, predict(lo, nx), col = 2, lwd = 2, lty = 2)
             do_poly(nx, pr$fit - cv * pr$se.fit, pr$fit + cv * pr$se.fit)
           }
           
           if (trend == 'gam') {
             k  <- if (missing(span))
               -1L else as.integer(span)
             gm <- mgcv::gam(y ~ s(x, k = -1L))
             pr <- mgcv::predict.gam(gm, data.frame(x = nx),
                                     type = 'response', se.fit = TRUE)
             lines(nx, pr$fit, col = 2, lwd = 2, lty = 2)
             do_poly(nx, pr$fit - cv * pr$se.fit, pr$fit + cv * pr$se.fit)
           }
         }
         
         panel.last
       }
  )
  
  title(xlab = xlab %||% sprintf('Average of %s and %s', xl, yl),
        ylab = ylab %||% sprintf('Difference between %s and %s', xl, yl))
  
  invisible(dd)
}
