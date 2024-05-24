hotelling.test = function(sample1, sample2, delta.0, alpha = 0.05, total_number_bonferroni_intervals = NULL){
  
  p = ncol(sample1)
  
  if(p != ncol(sample2)){
    stop("The number of variables in the two samples must be the same")
  }
  
  if(length(delta.0) != p){
    stop("The length of the vector delta.0 must be equal to the number of variables")
  }
  
  
  n1 = nrow(sample1)
  n2 = nrow(sample2)
  
  S1 = cov(sample1)
  S2 = cov(sample2)
  
  Sp = (n1 - 1)/(n1 + n2 - 2 ) * S1  + (n2 - 1)/ (n1 + n2 - 2) * S2
  
  x1 = sapply(sample1, mean)
  x2 = sapply(sample2, mean)
  
  
  T2 = t((x1 - x2 - delta.0))%*%solve((1/n1 + 1/n2) * Sp) %*% ( x1 - x2 - delta.0)
  cfr.fisher = qf(1-alpha, p, n1 + n2 - p - 1)
  
  
  
  P_val = 1 - pf(T2*(n1 + n2 - p - 1)/p/(n1+ n2 -2), p, n1 + n2 - p -1)
  
  if(is.null(total_number_bonferroni_intervals)){
    total_number_bonferroni_intervals = p
  }
  
  BonfCI = cbind(
    inf = (x1 - x2) - qt(alpha/2/total_number_bonferroni_intervals, n1 + n2 -2 )*sqrt((1/n1 + 1/n2) * diag(Sp)),
    center = (x1 - x2), 
    sup = (x1 - x2) + qt(alpha/2/total_number_bonferroni_intervals, n1 + n2 -2 )*sqrt((1/n1 + 1/n2) * diag(Sp))
  )
  
  
  print("Hotelling's T-squared test")
  cat("Statistic:", T2, "\n")
  cat("Reject H0:", T2 > cfr.fisher, "\n")
  cat("P-value:", P_val, "\n")
  cat("Method: Hotelling's T-squared test\n")
  cat("Bonferroni Confidence Intervals:\n")
  print(BonfCI)
  
  output <- list(
    statistic = T2,
    reject.H0 = T2 > cfr.fisher,
    p.value = P_val,
    BonfCI = BonfCI,
    method = "Hotelling's T-squared test"
  )
  class(output) <- "htest"
  
  return(output)
}