# log-linear models
# examples from von Eye book

melt <- function(data, nraters, weight = 'count') {
  tmp <- matrix(unlist(lapply(seq_len(nrow(data)),
                              function(x)
                                rep(data[x, 1:nraters],
                                    times = data[[weight]][x]))),
                ncol = nraters, byrow = TRUE)
return(tmp)
}

# pg 36
one <- data.frame(x = rep(1:4, each = 4),
                  y = 1:4,
                  count = c(4,3,1,0,1,5,2,2,
                            1,1,6,2,0,1,5,6))

one <- within(one, 
              delta <- ifelse(x == y, 1, 0))
one.m <- melt(one, 2)

# pg 65
two <- data.frame(r1 = rep(1:3, each = 9),
                  r2 = rep(1:3, each = 3),
                  r3 = 1:3,
                  count = c(4,3,6,2,1,3,2,2,17,
                            0,1,2,1,1,1,0,0,4,
                            0,1,3,0,1,8,0,4,96))

two <- within(two, {
  delta123 <- ifelse(r1 == r2 & r2 == r3, 1, 0)
  delta23 <- ifelse(r2 == r3, 1, 0)
  delta13 <- ifelse(r1 == r3, 1, 0)
  delta12 <- ifelse(r1 == r2, 1, 0)
})

two.m <- melt(two, 3)

three <- data.frame(a = rep(1:3, each = 3),
                    b = 1:3,
                    count = c(17,27,3,16,45,14,1,3,3))

three <- within(three, {
  trend <- ifelse(a == b, 0, ifelse(a > b, -1, 1))
  delta <- ifelse(a == b, 1, 0)
})

three.m <- melt(three, 2)



### kappa 
# library(psych)
# cohen.kappa(matrix(one$count, nrow = 4, byrow = TRUE))
library(irr)
kappa2(one.m)
kappam.fleiss(two.m)
# matching sas controlling output
kappa2(two.m[two.m[ , 1] == 1, 2:3])
kappa2(two.m[two.m[ , 1] == 2, 2:3])
kappa2(two.m[two.m[ , 1] == 3, 2:3])
kappa2(two.m[ , 2:3])
kappa2(three.m)
kappa2(three.m, weight = 'equal')



### log-linear models
?loglin
MASS:::loglm1.default
library(MASS)
# base model, main effects only
one1 <- loglm(count ~ x + y, data = one)
# equal weight agreement model
one2 <- loglm(count ~ x + y + delta, data = one)
resid(one2)
fitted(one2)

fit1 <- glm(count ~ factor(x) + factor(y), data = one, family = 'poisson')
fit2 <- glm(count ~ x + y + delta, data = one, family = 'poisson')


# pg 65, section 2.1
# base model, main effects only 
two1 <- loglm(count ~ r1 + r2 + r3, data = two)
# simultaneous agreement between pairs of raters (section 2.4.1.1)
two2 <- loglm(count ~ r1 + r2 + r3 + delta12 + delta13 + delta23, data = two)
# agreement among all raters (section 2.4.1.2)
two3 <- loglm(count ~ r1 + r2 + r3 + delta12 + delta13 + delta23 + delta123, 
              data = two)
# agreement among all raters, include only the delta for all raters 
# (sect 2.4.1.2)
two4 <- loglm(count ~ r1 + r2 + r3 + delta123, data = two)


# other models in section 2.4.2 for rater-specific trends
# base model, main effects only
three1 <- loglm(count ~ a + b, data = three)
# equal weight agreement model with trend (section 2.4.1.1)
three2 <- loglm(count ~ a + b + delta + trend, data = three)
three2[c('deviance','nobs','df','terms')]