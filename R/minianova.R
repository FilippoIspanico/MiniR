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

#' Multivariate Analysis of Variance Homogeneity and Normality Test
#'
#' This function performs homogeneity and normality tests for multivariate analysis of variance (MANOVA) 
#' or univariate analysis of variance (ANOVA) designs. It can handle one-way and two-way designs for 
#' both univariate and multivariate data.
#'
#' @param data A data frame or vector containing the response variable(s).
#' @param effect1 A vector specifying the levels of the first factor.
#' @param effect2 An optional vector specifying the levels of the second factor for two-way designs. Default is NULL.
#'
#' @return The function returns different objects based on the input:
#'   \itemize{
#'     \item For univariate data with one factor: A matrix with Shapiro-Wilk test p-values and variances for each level of effect1.
#'     \item For multivariate data with one factor: A vector of Mardia's multivariate normality test p-values for each level of effect1.
#'     \item For univariate data with two factors: A matrix with Shapiro-Wilk test p-values and variances for each combination of effect1 and effect2.
#'     \item For multivariate data with two factors: A vector of Mardia's multivariate normality test p-values for each combination of effect1 and effect2.
#'   }
#'
#' @details 
#' The function uses Shapiro-Wilk test for univariate normality and Mardia's test (from the MVN package) 
#' for multivariate normality. It also calculates variances (for univariate data) or covariance matrices 
#' (for multivariate data) for each group.
#' 
#' For multivariate data, the function produces plots of the covariance matrices.
#'
#' @note 
#' This function requires the 'MVN' package for multivariate normality tests. 
#' Make sure it's installed before using this function for multivariate data.
#'
#' @export
#' @export m.aov.hptest
#'
#' @examples
#' \dontrun{
#' # One-way univariate example
#' data <- rnorm(100)
#' effect1 <- rep(c("A", "B"), each = 50)
#' result1 <- m.aov.hptest(data, effect1)
#'
#' # Two-way multivariate example
#' data <- data.frame(y1 = rnorm(120), y2 = rnorm(120))
#' effect1 <- rep(c("X", "Y", "Z"), each = 40)
#' effect2 <- rep(rep(c("Low", "High"), each = 20), 3)
#' result2 <- m.aov.hptest(data, effect1, effect2)
#' }
#'
#' @importFrom graphics par image
#' @importFrom stats var cov shapiro.test
#' @importFrom MVN mvn

m.aov.hptest = function(data, effect1, effect2 = NULL){
  if (!requireNamespace("MVN", quietly = TRUE)) {
    stop("Package 'MVN' is required. Please install it.")
  }
  
  if(class(data) == "data.frame"){
    p = ncol(data)
  }else{
    data = data.frame(data)
    p = 1
  }
  
  
  treat1 = unique(effect1)
  treat2 = unique(effect2)
  
  p.values = NULL
  S = NULL
  
  
  if(is.null(effect2) & p == 1){
    
    for(i in 1:length(treat1)){
      indexes = which(effect1 == treat1[i])
      data.subset = data[indexes, ]
      p.values = cbind(
        p.values, 
        shapiro.test(data.subset)$p.value
      )
      S  = cbind(S, var(data.subset))
    }
    
    result = rbind(
      p.values, S)
    rownames(result) = c("shapiro p-value", "variance")
    
    return(result)
  }
  
  if(is.null(effect2)){
    par(mfrow = c(1, length(treat1)))
    for(i in 1:length(treat1)){
      indexes = which(effect1 == treat1[i])
      data.subset = data[indexes, ]
      
      p.values = cbind(
        p.values,
        MVN::mvn(data.subset)$multivariateNormality$`p value`
      )
      effect_label = as.character(treat1[i])
      cat("covariance matrix of effect: ",effect_label, "\n")
      image(cov(data.subset), main = paste("cov", effect_label))
      print(cov(data.subset))
      cat("\n")
    }
    par(mfrow = c(1, 1))
    rownames(p.values) = "mvn - p.value"
    colnames(p.values) = treat1
    return(p.values)
  }
  
  if(p == 1){
    for(i in 1:length(treat1)){
      for(j in 1:length(treat2)){
        
        indexes = which(effect1 == treat1[i] & effect2 == treat2[j])
        data.subset = data[indexes, ]
        p.values = cbind(
          p.values, 
          shapiro.test(data.subset)$p.value
        )
        S  = cbind(S, var(data.subset))
      }
    }
    
    result = rbind(
      p.values, S)
    rownames(result) = c("shaprio p-value", "variance")
    
    return(result)
  }
  
  par(mfrow = c(length(treat1), length(treat2)))
  for(i in 1:length(treat1)){
    for(j in 1:length(treat2)){
      
      indexes = which(effect1 == treat1[i] & effect2 == treat2[j])
      data.subset = data[indexes, ]
      
      p.values = cbind(
        p.values,
        MVN::mvn(data.subset)$multivariateNormality$`p value`
      )
      effect_label1 = as.character(treat1[i])
      effect_label2 = as.character(treat2[j])
      
      cat("covariance matrix of effects: ",effect_label1, "-", effect_label2, "\n")
      print(cov(data.subset))
      image(cov(data.subset), main = paste("cov", effect_label1, "-", effect_label2))
      cat("\n")
      
    }
  }
  par(mfrow = c(1, 1))
  rownames(p.values) = "mvn - p.value"
  print(interaction(treat1, treat2))
  colnames(p.values) = interaction(treat1, treat2)
  return(p.values)
}



