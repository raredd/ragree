# interrater agreement
# rmac package

# rmac method is insensitive to changes in marginal probabilities
# Fay says this is an advantage to his method

library(rmac)
?rmac
?wkappa

#The Multiple Sclerosis Diagnoses Example (from Fay 2005)
#the original data
msd1<- as.table(matrix(data = c(38,5,0,1,33,11,3,0,10,14,5,6,3,7,3,10), 4,4, byrow = TRUE))
msd1
#the data with cell counts (1,3) and (3,1) reversed
msd2<- as.table(matrix(data = c(38,5,10,1,33,11,3,0,0,14,5,6,3,7,3,10), 4,4, byrow = TRUE))
msd2

# fmac = c(.208,.186)
# rmac = c(.178,.178)



# with balanced marginals (as much as possible)
msd1<- as.table(matrix(data = c(38,19,5,2,19,11,9,4,5,8,5,5,2,3,4,10), 4,4, byrow = TRUE))
msd1
# with odd values reversed
msd2<- as.table(matrix(data = c(38,19,5,2,19,11,9,3,5,8,5,4,2,4,5,10), 4,4, byrow = TRUE))
msd2

# fmac = c(.178,.178)
# rmac = c(.178,.178)



# upper/lower triangular with same total
msd1<- as.table(matrix(data = c(38,38,10,4,0,11,17,7,0,0,5,9,0,0,0,10), 4,4, byrow = TRUE))
msd1
# with odd values reversed
msd2<- as.table(matrix(data = c(38,38,10,4,0,11,17,7,0,0,5,9,0,0,0,10), 4,4, byrow = FALSE))
msd2

# fmac = c(.224,.224)
# rmac = c(.178,.178)



# with balanced marginals
msd1<- as.table(matrix(data = c(38,38,7,7,0,11,19,5,0,0,5,9,0,0,0,10), 4,4, byrow = TRUE))
msd1
# with odd values reversed
msd2<- as.table(matrix(data = c(38,39,6,7,0,11,19,5,0,0,5,9,0,0,0,10), 4,4, byrow = TRUE))
msd2

# fmac = c(.224,.223)
# rmac = c(.178,.177)





#calculate the FMAC of each data set
wkappa(msd1, method = "fmac")
wkappa(msd2, method = "fmac")

#calculate the FMAC of each data set
wkappa(msd1, method = "rmac")
wkappa(msd2, method = "rmac")

#