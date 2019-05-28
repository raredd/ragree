### misc
# print.ragree, print.ragree1, lr_test
###


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
    out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 4))))
  if (!is.null(x$parameter)) 
    out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 3))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = digits)
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<")
      fp else paste("=", fp)))
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
#' lr_test(dat)
#' 
#' ## compare to
#' loglin(dat, margin = 1:2)
#' 
#' @export

lr_test <- function(x, print = TRUE) {
  chi.out <- chisq.test(x, correct = FALSE)
  dat     <- chi.out$observed
  ratios  <- dat / chi.out$expected
  
  ## lr chisq statistic, degrees of freedom, p-value
  lrt  <- 2 * sum(dat * log(ratios))
  df   <- unname(chi.out$parameter)
  pval <- pchisq(lrt, df, lower.tail = FALSE)
  
  if (print) 
    print(sprintf('LR: %s on %s degrees of freedom (p-value: %s)',
                  round(lrt, digits = 2L), df, pvalr(pval)))
  
  list(lrt = lrt, df = df, p.value = pval, pearson = unname(chi.out$statistic))
}
