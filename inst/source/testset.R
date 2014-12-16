# testset.sas
# magree.mac


library(irr)
library(plyr)
library(reshape2)


dat <- read.csv('/users/rawr/documents/dfci/edie/rater agreement/testset.csv',
                header = TRUE, stringsAsFactors = FALSE)

## dat has one "marker"
kappam.fleiss(dcast(dat, patient ~ rater)[ , -1], exact = FALSE)$value

## dat2 has multiple "markers"
# simulate kappa data
kap <- function(x, y, z) {
  kap <- melt(replicate(108 / 4, 
                        c(sample(c('mild','mod','sev'), 4, 
                                 replace = TRUE, prob = c(x, y, z)))))$value
  return(as.character(kap))
}

set.seed(1618)
dat2 <- within(dat, {
  marker4 <- kap(0,.1,.9)
  marker3 <- kap(0,.5,.5)
  marker2 <- kap(.3,.3,.3)
})
names(dat2)[3] <- 'marker1'

# long way
kappam.fleiss(dcast(dat2[-(4:6)], patient ~ rater)[-1], exact = FALSE)$value
kappam.fleiss(dcast(dat2[-c(3,5:6)], patient ~ rater)[-1], exact = FALSE)$value
kappam.fleiss(dcast(dat2[-c(3:4,6)], patient ~ rater)[-1], exact = FALSE)$value
kappam.fleiss(dcast(dat2[-c(3:5)], patient ~ rater)[-1], exact = FALSE)$value


(dcast(tmp, patient ~ rater)[-1])

dat2.s <- split(dat2, dat2$rater)

markers <- paste0('marker',1:4)

sapply(markers, function(x)
  kappam.fleiss(dcast(dat2[c('patient', 'rater', x)],
                      patient ~ rater, value.var = x)[ , -1], 
                exact = FALSE)$value)



# testset.sas

xtabs( ~ patient + category, data = dat)
with(dat, table(category) / sum(table(category)))

# patval
patval <- ddply(dat, .(patient), summarise, 
                count = length(patient),
                percent = length(patient) / nrow(dat))

# outval2
outval2 <- ddply(dat, .(category), summarise,
                 count = length(category),
                 percent = length(category) / nrow(dat) * 100)

# outval1
outval1 <- unique(ddply(dat[-2], .(patient, category), transform, 
                       njk = length(category),
                       njk2 = length(category)**2))

# outval
outval <- merge(patval, outval1, by = 'patient')
outval <- within(outval, tempnk <- njk * (count - 1))

# final 
final <- ddply(outval, .(category), summarise,
              freq = length(category),
              sumnjk = sum(njk),
              sumnjk2 = sum(njk2),
              poss = sum(tempnk))

# outval2
outval2 <- within(outval2, {
  pj <- percent / 100
  pj2 <- (percent / 100) ** 2
  pj3 <- (percent / 100) ** 3
})

# final
final <- merge(final, outval2, by = 'category')
final <- within(final, {
  actualj <- sumnjk2 - sumnjk
  pej <- pj2
  agreej <- actualj / poss
  kappaj <- (agreej - pj) / (1 - pj)
})

# meanval
meanval <- with(final, data.frame(actual = sum(actualj),
                                  pe = sum(pj2),
                                  p3 = sum(pj3),
                                  poss = sum(poss)))

# meanval
meanval <- within(meanval, {
  agree <- actual / poss
  kappa <- ifelse(pe != 1, (agree - pe) / (1 - pe), 1)
  nr <- 4
  c1 <- 2 * nr - 3
  c2 <- 2 * (nr - 2)
  pe2 <- pe ** 2
  # variance of generalized kappa (page 279 woolson)
  vk <- ifelse(pe != 1, 
               (2 * (pe - c1 * pe2 + c2 * p3)) / (poss * (1 - pe) ** 2), 
               0)
  se <- sqrt(vk)
})