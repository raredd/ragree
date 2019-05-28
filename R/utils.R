### utils
# %||%, pvalr, roundr
###


`%||%` <- function(x, y) {
  if (is.null(x) || !nzchar(x))
    y else x
}

pvalr <- function(pv, sig.limit = 0.001, digits = 3L, html = FALSE,
                  show.p = FALSE, journal = TRUE, ...) {
  ## rawr::pvalr
  stopifnot(
    sig.limit > 0,
    sig.limit < 1
  )
  
  show.p <- show.p + 1L
  html   <- html + 1L
  
  sapply(pv, function(x, sig.limit) {
    if (is.na(x) | !nzchar(x))
      return(NA)
    if (x >= 0.99)
      return(paste0(c('', 'p ')[show.p], c('> ', '&gt; ')[html], '0.99'))
    if (x >= 0.9 && !journal)
      return(paste0(c('', 'p ')[show.p], c('> ', '&gt; ')[html], '0.9'))
    if (x < sig.limit) {
      paste0(c('', 'p ')[show.p], c('< ', '&lt; ')[html],
             format.pval(sig.limit, ...))
    } else {
      nd <- if (journal)
        c(digits, 2L)[findInterval(x, c(-Inf, .1, Inf))]
      else c(digits, 2L, 1L)[findInterval(x, c(-Inf, .1, .5, Inf))]
      paste0(c('', 'p = ')[show.p], roundr(x, nd))
    }
  }, sig.limit)
}

roundr <- function(x, digits = 1) {
  ## rawr::roundr
  res <- sprintf(paste0('%.', digits, 'f'), x)
  zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
  res[res == paste0('-', zzz)] <- zzz
  res
}
