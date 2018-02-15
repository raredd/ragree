context('unalike')

set.seed(1)
x <- sample(letters, 1000L, replace = TRUE)

test_that('unalike1: table and vector inputs are interchangeable', {
  expect_identical(
    ragree:::unalike1(x),
    ragree:::unalike1(table(x))
  )
  
  expect_identical(
    unalike(x, method = 1L),
    unalike(table(x), method = 1L)
  )
})


test_that('unalike2: table and vector inputs are interchangeable', {
  expect_identical(
    ragree:::unalike2(x),
    ragree:::unalike2(table(x))
  )
  
  expect_identical(
    unalike(x, method = 2L),
    unalike(table(x), method = 2L)
  )
})

set.seed(1)
mat <- replicate(10, sample(letters, 1000L, replace = TRUE))
colnames(mat) <- paste('rater', 1:10)

dat <- data.frame(
  id    = rep(seq.int(nrow(mat)), ncol(mat)),
  rater = rep(colnames(mat), each = nrow(mat)),
  score = c(mat)
)

test_that('unalike S3: matrix and data frame inputs are interchangeable', {
  expect_identical(
    unalike(mat)[1:7],
    unalike(dat)[1:7]
  )
})


rm(x, mat, dat)
