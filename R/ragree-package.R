#' Rater agreement
#' 
#' Various measures for quantifying rater agreement
#' 
#' @name ragree-package
#' @aliases ragree
#' @docType package
#' @title Rater agreement
#' @author \email{rredd@@jimmy.harvard.edu}
#' @keywords rater agreement kappa unalikeability unalike
NULL

#' Disease severity rating of 27 patients by four ophthalmologists per person
#' 
#' Sample data in which 27 patients with a particular retinopathy had the
#' severity of their retinopathy rated by four ophthalmologists.
#' 
#' @name disease.l
#' @docType data
#' @keywords data disease
#' @usage disease.l
#' @seealso \code{\link{disease.w}}
#' @format An object of class \code{data.frame} containing 108 observations and
#' 3 variables:
#' 
#' \tabular{rll}{
#' \tab \code{patient} \tab patient number \cr
#' \tab \code{rater} \tab rater number \cr
#' \tab \code{category} \tab disease severity \cr
#' }
#' 
NULL

#' Disease severity rating of 27 patients by four ophthalmologists per person
#' 
#' Sample data in which 27 patients with a particular retinopathy had the
#' severity of their retinopathy rated by four ophthalmologists.
#' 
#' @name disease.w
#' @docType data
#' @keywords data disease
#' @usage disease.w
#' @seealso \code{\link{disease.l}}
#' @format An object of class \code{data.frame} containing 27 observations and
#' 4 variables:
#' 
#' \tabular{rll}{
#' \tab \code{rater1} \tab opthalmologist 1 ratings \cr
#' \tab \code{rater2} \tab opthalmologist 2 ratings \cr
#' \tab \code{rater3} \tab opthalmologist 3 ratings \cr
#' \tab \code{rater4} \tab opthalmologist 4 ratings \cr
#' }
#' 
NULL

#' ragree class print method
#' @usage print(x, ...)
#' @param x object of class 'ragree'
#' @param ... other arguments passed to or from other methods
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