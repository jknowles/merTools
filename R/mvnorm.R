#' @title Fast rcpp rmvnorm sampler
#' @name rccp_rmvnorm
#' @param n number of samples
#' @param mu vector of means
#' @param sigma sigma matrix
#' @export
rcpp_rmvnorm <- function(n, mu, sigma)
{
  chol_S <- rcppeigen_get_chol(sigma)
  m <- dim(chol_S)[1]
  mat <- matrix(rnorm(n*m),nrow=m,ncol=n)
  if(class(mu) != "numeric"){
    mu <- as.numeric(mu)
  }
  return (t(chol_S%*%mat+mu))
}
