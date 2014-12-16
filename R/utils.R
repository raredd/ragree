#' ragree class print method
#' 
#' @param x object of class 'ragree'
#' @param ... other arguments passed to or from other methods
#' 
#' @export

print.ragree <- function (x, ...) {
  
  cat(' ', x$method, '\n\n', sep = '')
  cat(paste('  Subjects =', x$subjects, '\n'))
  cat(paste('    Raters =', x$raters, '\n'))
  cat(paste('Categories =', x$categories, '\n'))
  
  results <- paste(formatC(x$ragree.name, width = 10, flag = '+'), 
                   '=', format(x$value, digits = 3), '\n')
  cat(results)
  if (!is.null(x$statistic)) {
    statistic <- paste(formatC(x$stat.name, width = 10, flag = '+'), 
                       '=', format(x$statistic, digits = 3), '\n')
    cat('\n', statistic, sep = '')
    cat(paste('   p-value', ifelse(x$p.value < .001, 
                                   format.pval(x$p.value, eps = .001),
                                   paste0('= ', round(x$p.value, 3))),'\n'))
  }
  
  #   if (!is.null(x$detail)) {
  #     cat('\n')
  #     if (is.table(x$detail)) 
  #       print(x$detail)
  #     else print(x$detail, '\n')
  #   }
  
  if (!is.null(x$variance)) 
    cat('\n ', x$error, '\n', sep = '')
}

print.ragree1 <- function (x, digits = 4L, quote = TRUE, prefix = "", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$statistic)) 
    out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 
                                                              4))))
  if (!is.null(x$parameter)) 
    out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 
                                                              3))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = digits)
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == 
                                         "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to", 
                           less = "less than", greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char, 
            " ", x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n", 
        " ", paste(format(c(x$conf.int[1L], x$conf.int[2L])), 
                   collapse = " "), "\n", sep = "")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, ...)
  }
  cat("\n")
  invisible(x)
}

#' LR chi-squared test for two-way tables
#' 
#' Likelihood ratio chi-squared test for two-way contingency tables.
#' 
#' @param x two-way table
#' @param print logical; if \code{TRUE}, prints summary of test
#' @references 
#' \url{http://ww2.coastal.edu/kingw/statistics/R-tutorials/loglin.html}
#' 
#' @examples
#' data(Titanic)
#' dat <- margin.table(Titanic, c(2, 4))
#' lr.test(dat)
#' 
#' ## compare to
#' loglin(dat, margin = 1:2)
#' 
#' @export

lr.test <- function(x, print = TRUE) {
  chi.out <- chisq.test(x, correct = FALSE)
  dat <- chi.out[[6]]
  ratios <- dat / chi.out[[7]]
  
  ## lr chisq statistic, degrees of freedom, p-value
  lrt <- 2 * sum(dat * log(ratios))
  df <- unname(chi.out[[2]])
  pval <- pchisq(lrt, df, lower.tail = FALSE)
  
  if (print) 
    print(sprintf('LR: %s on %s degrees of freedom (p-value: %s)',
                  round(lrt, digits = 2), df, pvalr(pval)))
  
  return(list(lrt = lrt,
              df = df,
              p.value = pval,
              pearson = unname(chi.out[[1]])))
}

## rawr::pvalr
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE, 
         show.p = FALSE) {
  ## rawr::roundr
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  sapply(pvals, function(x, sig.limit) {
    if (is.na(x))
      return(NA)
    if (x >= 1)
      return('1')
    if (x < sig.limit) {
      if (show.p) p <- 'p ' else p <- ''
      if (html)
        return(sprintf('%s&lt; %s', p, format(sig.limit))) else
          return(sprintf('%s< %s', p, format(sig.limit)))
    } else {
      if (show.p) p <- 'p = ' else p <- ''
      if (x > .1)
        return(sprintf('%s%s', p, roundr(x, digits = 2))) else
          return(sprintf('%s%s', p, roundr(x, digits = digits)))
    }
  }, sig.limit = sig.limit)
}