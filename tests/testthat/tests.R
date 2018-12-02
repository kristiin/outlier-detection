test_that("upper interval is correct", {
  expect_equal(interval.upper(10,2,2.5),15)
})

test_that("lower interval is correct", {
  expect_equal(interval.lower(10,2,2.5),5)
})

#test_that("mad is calculated per each column"){
#  x<-madPerColumnFun(getTestData())
#  expect_equal(x,c(6,))
#}

getTestData <- function(){
  return (data.frame("VEL" = c(3,6,6,8,10), "HEIGHT" = c(11,43,22,21,15), "COHER" = c(0.6,0.3,0.8,0.9,0.7), stringsAsFactors = FALSE))

}
