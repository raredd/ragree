
(table.8.3 <- data.frame(expand.grid(
  a = c("Yes","No"),
  c = c("Yes","No"),
  m = c("Yes","No")),
  count = c(911,538,44,456,3,43,2,279)))

(dat <- within(table.8.3, {
  acm <- ifelse(a == c & c == m, 1, 0)
  cm <- ifelse(c == m, 1, 0)
  am <- ifelse(a == m, 1, 0)
  ac <- ifelse(a == c, 1, 0)
}))


library(MASS)
fitACM <- loglm(count ~ a * c * m, data = dat, param = TRUE, fit = TRUE) # ACM
fit1 <- glm(count ~ a * c * m, data = dat, family = 'poisson')
fitAC.AM.CM <- update(fitACM, . ~ . - a:c:m) # AC, AM, CM
fitAM.CM <- update(fitAC.AM.CM, . ~ . - a:c) # AM, CM
fitAC.M <- update(fitAC.AM.CM, . ~ . - a:m - c:m) # AC, M
fitA.C.M <- update(fitAC.M, . ~ . - a:c) # A, C, M


data.frame(table.8.3[-4],
           ACM = c(aperm(fitted(fitACM))),
           AC.AM.CM = c(aperm(fitted(fitAC.AM.CM))), 
           AM.CM = c(aperm(fitted(fitAM.CM))),
           AC.M = c(aperm(fitted(fitAC.M))), 
           A.C.M = c(aperm(fitted(fitA.C.M))))
