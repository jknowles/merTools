library(testthat)
library(merTools)

test_check("myPkg", filter = "^[a-m]")
test_check("myPkg", filter = "^[n-z]")
