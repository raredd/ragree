identity.weights <- function(categ) diag(length(categ))

quadratic.weights <- function(categ)
	outer(seq_along(categ), seq_along(categ), 
        function(x, y) 1 - (x - y) ** 2 / diff(range(x, y)) ** 2)

linear.weights <- function(categ) 
  outer(seq_along(categ), seq_along(categ), 
        function(x, y) 1 - abs(x - y) / diff(range(x, y)))

radical.weights <- function(categ) 
  outer(seq_along(categ), seq_along(categ), 
        function(x, y) 1 - sqrt(abs(x - y)) / sqrt(abs(diff(range(x, y)))))

ratio.weights <- function(categ)
  outer(seq_along(categ), seq_along(categ),
        function(x, y) 1 - ((x - y) / (x + y)) ** 2 / 
          (diff(range(categ)) / sum(range(categ))) ** 2)

circular.weights <- function(categ) {
  zzz <- outer(seq_along(categ), seq_along(categ),
               function(x, y) sin(pi * (x - y) / (diff(range(categ)) + 1)) ** 2)
  1 - zzz / max(zzz)
}

bipolar.weights <- function(categ) {
  zzz <- outer(seq_along(categ), seq_along(categ),
               function(x, y)
                 ifelse(x == y, 0, (x - y) ** 2 / (((x + y) - 2 * min(x)) * 
                          (2 * max(x) - (x + y)))))
  1 - zzz / max(zzz)
}

ordinal.weights <- function(categ) {
  zzz <- outer(seq_along(categ), seq_along(categ),
               function(x, y) {
                 nkl <- pmax(x, y) - pmin(x, y) + 1
                 nkl * (nkl - 1) / 2})
  1 - zzz / max(zzz)
}
