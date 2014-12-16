## data sets in ragree package

## long format
disease.l <- read.csv('./source/testset.csv',
                      header = TRUE, stringsAsFactors = FALSE)

## wide format using ratings argument
disease.w <- reshape(disease.l, timevar = 'rater', idvar = 'patient',
                     direction = 'wide')[-1]
names(disease.w) <- c('rater1','rater2','rater3','rater4')


save(disease.l, file = './data/disease.l.rda')
save(disease.w, file = './data/disease.w.rda')