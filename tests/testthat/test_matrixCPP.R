###############################################
context("Test matrix dimensions")
################################################
#
# tests <- expand.grid(n=c(2, 10, 100, 1000), dim=c(1, 2, 10, 100))
# set.seed(12451)
# for(i in 1:nrow(tests)){
#   mu <- rnorm(tests[i, "dim"], 0, 1)
#   w <- runif(tests[i, "dim"])
#   k <- length(w)
#   if(length(w) == 1){
#     w <-w
#   } else{
#     w <- diag(w)
#   }
#   x <- matrix(rnorm(k^2), nrow=k, ncol=k) %*% w
#   x <- x/sqrt(rowSums(x^2))
#   a <- x %*% t(x)
#   p <- FastGP::rcpp_rmvnorm_stable(n = tests[i, "n"], mu = mu, S = a)
#   if(tests[i, "dim"] == 1){
#     expect_is(p, "numeric")
#     expect_equal(length(p), tests[i, "n"])
#   } else{
#     expect_is(p, "matrix")
#     expect_equal(nrow(p), tests[i, "n"])
#     expect_equal(ncol(p), tests[i, "dim"])
#   }
# }


# library(microbenchmark)
# #
# microbenchmark(
#   for(i in 1:nrow(tests)){
#     mu <- rnorm(tests[i, "dim"], 0, 1)
#     w <- runif(tests[i, "dim"])
#     k <- length(w)
#     if(length(w) == 1){
#       w <-w
#     } else{
#       w <- diag(w)
#     }
#     x <- matrix(rnorm(k^2), nrow=k, ncol=k) %*% w
#     x <- x/sqrt(rowSums(x^2))
#     a <- x %*% t(x)
#     # p <- mvrnormArma(tests[i, "n"], mu = mu, sigma = a)
#     p <- mvtnorm::rmvnorm(tests[i, "n"], mean =  mu, sigma = a)
#     # p <- FastGP::rcpp_rmvnorm(n = tests[i, "n"], mu = mu, S = a)
#   }, times = 50
# )


