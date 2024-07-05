#' Calculate Confidence Interval for Effect Difference
#'
#' This function calculates the confidence interval for the difference between two effects
#' using the t-distribution and Bonferroni correction.
#'
#' @param effect1 Integer. Index of the first effect.
#' @param effect2 Integer. Index of the second effect.
#' @param alpha Numeric. Significance level (e.g. 0.05 for 95 percent confidence).
#' @param k Integer. Total number of intervals being calculated (for Bonferroni correction).
#' @param g Integer. Number of groups or treatments.
#' @param n Numeric vector. Sample sizes for each group.
#' @param mu Numeric vector. Means of each group.
#' @param Spool Numeric. Pooled variance estimate.
#' @param effect.names Character vector. Names of the effects or treatments.
#'
#' @return A matrix with one row, containing the lower (inf) and upper (sup) bounds of the confidence interval.
#'         The row is named as the difference between the two effects.
#'
#' @examples
#' ci.diff.effect(1, 2, 0.05, 3, 3, c(10, 10, 10), c(5, 6, 7), 2, c("A", "B", "C"))
#' @export ci.diff.effect
#' @export
ci.diff.effect = function(effect1, effect2, alpha, k, g, n, mu, Spool, effect.names){
  
  quantile = qt(1-alpha/2/k, sum(n)-g)
  correction = sqrt( (1/n[effect1] + 1/n[effect2]) * Spool )
  
  CI = cbind(
    inf = mu[effect1] - mu[effect2] - quantile*correction, 
    sup = mu[effect1] - mu[effect2] + quantile*correction
  )
  
  rownames(CI) = paste(effect.names[effect1], "-", effect.names[effect2])
  return(CI)
}


#' Calculate ANOVA Bonferroni Confidence Intervals for Differences
#'
#' This function calculates Bonferroni-corrected confidence intervals for all pairwise
#' differences between treatment effects in an ANOVA model.
#'
#' @param fit.aov An object of class 'aov', typically the result of a call to aov().
#' @param effect.factor Factor. The factor variable representing all treatments.
#' @param mu Numeric vector. Means of each treatment group.
#' @param alpha Numeric. Significance level (e.g., 0.05 for 95% confidence).
#' @param k Integer, optional. Number of intervals to calculate. If NULL (default),
#'        all pairwise differences are calculated.
#' @param round.3 Logical. If TRUE (default), round the results to 3 decimal places.
#'
#' @return A matrix where each row represents a pairwise difference between treatments.
#'         The columns give the lower (inf) and upper (sup) bounds of each confidence interval.
#'
#' @examples
#' # Assuming 'my_aov' is an aov object and 'treatment' is your factor variable
#' aov.bci.diff(my_aov, treatment, c(10, 12, 15), 0.05)
#'
#' @note If k is manually set, it should only be done when the problem specifically
#'       asks to compute a subset of all differences. For all differences, leave k as NULL.
#' @export aov.bci.diff
#' @export
aov.bci.diff  = function(fit.aov, effect.factor, mu, alpha, k = NULL, round.3 = TRUE){
  
  stopifnot("Input effect.factor must be the factor of all treatments. It must be of class factor. " = is.factor(effect.factor))
  
  n = table(effect.factor)
  treat = levels(effect.factor) 
  g = length(treat)
  
  if( !is.null(k) )
  {
    warning("you manually set the number of intervals k!
this is correct only if the text of the problem ask you to compute a subset of all 
the differnces.
If you want to compute all differnces remove the input k from the function call!")
    
  }else{
    k = g*(g-1)/2
  }
  
  S = sum(fit.aov$residuals^2)/fit.aov$df.residual
  CI = NULL
  
  for(i in 1:(g-1))
    for(j in (i+1):g){
      CI = rbind(
        CI, 
        ci.diff.effect(i, j, alpha, k, g, n, mu, S, effect.names = treat)
      )
    }
  
  if(round.3)
    CI = round(CI, digits = 3)
  
  
  return(CI)
}



