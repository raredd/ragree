### utils
# %||%, roundr, pvalr, recycle
###

'%||%' <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

roundr <- function(x, digits = 1) {
  ## rawr::roundr
  res <- sprintf(paste0('%.', digits, 'f'), x)
  zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
  res[res == paste0('-', zzz)] <- zzz
  res
}

pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE,
                  show.p = FALSE) {
  ## rawr::pvalr
  sapply(pvals, function(x, sig.limit) {
    if (is.na(x)) return(NA)
    if (x >= 1) return('> 0.99')
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

recycle <- function(x, y) {
  ## iplotr:::recycle
  ## for a vector y, repeats y to length of x
  ## recycle(1:5, 1:2)
  ## recycle(1:2, c('red','blue','green'))
  lx <- length(x)
  ly <- length(y)
  rep(y, ceiling(lx / ly))[seq_along(x)]
}
