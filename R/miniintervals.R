#' Simultaneous Confidence Intervals based on Hotelling's T-squared distribution
#'
#' This function calculates the simultaneous confidence intervals (SCI) for linear combinations of means 
#' based on Hotelling's T-squared distribution.
#'
#' @param x.mean A numeric vector representing the sample means.
#' @param S A numeric matrix representing the sample covariance matrix.
#' @param alpha A numeric value representing the significance level (e.g., 0.05 for 95% confidence).
#' @param A A numeric matrix where each row represents a linear combination of the means.
#'
#' @return A numeric matrix with three columns: `inf` for the lower bound of the SCI, `center` for the 
#' point estimates of the linear combinations, and `sup` for the upper bound of the SCI. Each row corresponds 
#' to a linear combination defined by a row of the matrix `A`.
#'
#' @details This function computes the simultaneous confidence intervals using Hotelling's T-squared distribution. 
#' The critical value is obtained from the Fisher's F distribution. The simultaneous confidence intervals are 
#' computed for the linear combinations specified by the rows of the matrix `A`.
#'
#' @examples
#' x.mean <- c(5, 10)
#' S <- matrix(c(4, 2, 2, 3), nrow = 2)
#' alpha <- 0.05
#' A <- matrix(c(1, 0, 0, 1), nrow = 2)
#' sim.T2(x.mean, S, alpha, A)
#'
#' @importFrom stats qf
#' @export

sim.T2 = function(x.mean, S, alpha, A){
  
  cfr.fisher = (n-1)*p/(n-p)*qf(1-alpha, p, n-p)
  
  SCI = cbind(
    inf = A%*%x.mean - sqrt(cfr.fisher) * sqrt(diag(A%*%S%*%t(A))/n),
    center = A%*%x.mean,
    sup = A%*%x.mean + sqrt(cfr.fisher) * sqrt(diag(A%*%S%*%t(A))/n)
  )
  
  return(SCI)
  
}



#' Compute Bonferroni-T2 Confidence Intervals
#'
#' This function computes Bonferroni-T2 confidence intervals for linear combinations of means.
#'
#' @param x.mean A numeric vector representing the sample means.
#' @param S A numeric matrix representing the sample covariance matrix.
#' @param alpha A numeric value representing the significance level for the confidence intervals.
#' @param A A numeric matrix where each row represents a linear combination of the means for which the confidence interval is computed.
#' @param k An optional integer representing the number of comparisons. If `NULL`, it defaults to the number of rows of `A`.
#'
#' @details 
#' The function calculates the Bonferroni-adjusted confidence intervals for the linear combinations of the means specified in the rows of `A`. The adjustment accounts for multiple comparisons to control the overall type I error rate.
#' 
#' If the provided `k` is greater than the number of rows in `A`, a warning is issued suggesting manual adjustments for other Bonferroni intervals. If `k` is smaller than the number of rows in `A`, a warning is issued and `k` is set to the number of rows in `A`.
#'
#' @return A matrix with three columns:
#' \describe{
#'   \item{inf}{The lower bounds of the confidence intervals.}
#'   \item{center}{The point estimates for the linear combinations of the means.}
#'   \item{sup}{The upper bounds of the confidence intervals.}
#' }
#'
#' @importFrom stats qt
#' @examples
#' x.mean <- c(2.5, 3.0, 4.5)
#' S <- matrix(c(1, 0.5, 0.3,
#'               0.5, 2, 0.4,
#'               0.3, 0.4, 3), nrow = 3)
#' alpha <- 0.05
#' A <- matrix(c(1, 0, 0,
#'               0, 1, 0,
#'               0, 0, 1), nrow = 3)
#' bonf.T2(x.mean, S, alpha, A)
#'
#' @export


bonf.T2 = function(x.mean, S, alpha, A, k = NULL){
  if(is.null(k))
    k  = nrow(A)
  else
  {
    if( k > nrow(A) )
      warning(paste("Provided k (", k, ") is greater than nrow(A) (", nrow(A), "). 
                    This is correct iff you are computing other Bonferroni intervals manually!", sep=""))
    if( k < nrow(A))
      warning(paste("Provided k (", k, ") is smaller than nrow(A) (", nrow(A), "). 
                    This is not correct! Proceeding with computing k = nrow(A) (", nrow(A), 
                    ") Bonferroni confidence intervals. If this is not your intention, please modify matrix A!", sep=""))
  }
  
  cfr.t = qt(1-alpha/2/k, n - 1)
  
  BCI = cbind(
    inf = A%*%x.mean - cfr.t * sqrt(diag(A%*%S%*%t(A))/n), 
    center = A%*%x.mean,
    sup = A%*%x.mean + cfr.t * sqrt(diag(A%*%S%*%t(A))/n)
  )
  
  return(BCI)
  
}