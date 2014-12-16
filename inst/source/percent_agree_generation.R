## data generation for desired percent agreement

set.seed(1618)

pct_agree <- .95
x <- sample(c('N', 'P'), 500, replace = TRUE)

sampleneg <- function() 
  sample(c('N', 'P'), 1, prob = c(pct_agree, 1 - pct_agree))
samplepos <- function() 
  sample(c('N', 'P'), 1, prob = c(1 - pct_agree, pct_agree))

y <- sapply(1:length(x),
            function(ii) ifelse(x[ii] == 'P', samplepos(), sampleneg()))

z <- sapply(1:length(x),
            function(ii) ifelse(x[ii] == 'P', samplepos(), sampleneg()))

cd10 <- data.frame(a = x, b = y, c = z, stringsAsFactors = FALSE)
(c_table <- table(cd10[1:3]))
sum(cd10[1] == cd10[2]) / nrow(cd10)
sum(cd10[1] == cd10[3]) / nrow(cd10)


c_table <- cbind(expand.grid(c('N','P'), c('N','P'), c('N','P')), c(c_table))
names(c_table) <- c('a','b','c','count')

c_table <- within(c_table, {
  abc <- ifelse(a == b & b == c, 1, 0)
  bc <- ifelse(b == c, 1, 0)
  ac <- ifelse(a == c, 1, 0)
  ab <- ifelse(a == b, 1, 0)
})

## independent model
fit1 <- glm(count ~ factor(a) + factor(b) + factor(c),
           data = c_table, family = 'poisson')
summary(fit1)

## independent model + pairwise
fit2 <- glm(count ~ factor(a) + factor(b) + factor(c) + ab + ac + bc,
           data = c_table, family = 'poisson')
summary(fit2)

## independent model + pairwise + three-way
fit3 <- glm(count ~ factor(a) + factor(b) + factor(c) + 
             ab + ac + bc + abc,
           data = c_table, family = 'poisson')
summary(fit3)
## Coefficients: (1 not defined because of singularities)