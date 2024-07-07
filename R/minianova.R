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
  column_names = character(0)  
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
      
      column_names = c(column_names, paste(effect_label1, effect_label2, sep = "-"))
    }
  }
  par(mfrow = c(1, 1))
  rownames(p.values) = "mvn - p.value"

  colnames(p.values) = column_names
  return(p.values)
}

#' Custom version of colMeans
#'
#' colMeans but works also when p=1
#'
#' @param data A data frame or vector .
#' 
#' @export
#' @export colMeans
#'

colMeans = function(data){
  
  p = ncol(data)
  if(is.null(p))
  {
    result = matrix(mean(data), nrow = 1, ncol = 1)
    return(result)
  }
  
  if(p == 1){
    result = matrix(mean(data[,1]), nrow = 1, ncol = 1)
    return(result)
  }
  
  return(base::colMeans(data))
}



tau = function(data, treat, effect, treat2 = NULL, effect2 = NULL){
  
  if(is.null(treat2) | is.null(effect2))
    return(MiniR::colMeans(data[which(effect == treat) ,]))
  return(MiniR::colMeans(data[which(effect == treat & effect2 == treat2) ,]))
  
}


CI.diff = function(tau1, tau2, pos, n1, n2, k.bonf, alpha, fit.manova, p){
  
  if(p != 1){
    Spool = summary(fit.manova)$SS$Residuals
    Spool = Spool/fit.manova$df.residual
  }else{
    Spool = matrix(0, nrow = 1, ncol = 1)
    Spool[1, 1] = sum(fit.manova$residuals^2)/fit.manova$df.residual
  }
  
  quantile = qt(1-alpha/2/k.bonf, fit.manova$df.residual)
  correction = sqrt(Spool[pos, pos]*(1/n1 + 1/n2))
  
  ci = cbind(
    inf = tau1[pos] - tau2[pos] - quantile*correction,
    sup = tau1[pos] - tau2[pos] + quantile*correction
  )
  
  return(ci)
  
}
# n1 & n2 are not the observations per group! 
# they are the number of observations that are used to compute tau1 & tau2!



CI.effect = function(tau1, pos, n1, x.mean, n.tot, k.bonf, alpha, fit.manova, p){
  
  if(p != 1){
    Spool = summary(fit.manova)$SS$Residuals
    Spool = Spool/fit.manova$df.residual
  }else{
    Spool = matrix(0, nrow = 1, ncol = 1)
    Spool[1, 1] = sum(fit.manova$residuals^2)/fit.manova$df.residual
  }
  
  quantile = qt(1-alpha/2/k.bonf, fit.manova$df.residual)
  correction = sqrt(Spool[pos, pos]*(1/n1 + 1/n.tot))
  
  ci = cbind(
    inf = tau1[pos] - x.mean[pos] - quantile*correction,
    sup = tau1[pos] - x.mean[pos] + quantile*correction
  )
  
  return(ci)
  
}



CI = function(tau1, pos, n1, k.bonf, alpha, fit.manova, p){
  
  if(p != 1){
    Spool = summary(fit.manova)$SS$Residuals
    Spool = Spool/fit.manova$df.residual
  }else{
    Spool = matrix(0, nrow = 1, ncol = 1)
    Spool[1, 1] = sum(fit.manova$residuals^2)/fit.manova$df.residual
  }
  
  quantile = qt(1-alpha/2/k.bonf, fit.manova$df.residual)
  correction = sqrt(Spool[pos, pos]*(1/n1))
  
  ci = cbind(
    inf = tau1[pos] - quantile*correction,
    sup = tau1[pos] + quantile*correction
  )
  
  return(ci)
  
  
  
}




#' Compute Bonferroni Intervals for MANOVA or ANOVA
#'
#' This function calculates Bonferroni intervals for Multivariate Analysis of Variance (MANOVA) 
#' or univariate Analysis of Variance (ANOVA). It can handle one-way or two-way designs.
#' The function provides confidence intervals for differences between treatment effects 
#' and for the effects themselves.
#'
#' @param fit.manova An object representing a fitted MANOVA or ANOVA model.
#' @param data A dataframe containing the response variable(s).
#' @param effects1 A factor representing the first treatment effect.
#' @param effect1.name A character string for the name of the first effect (currently unused in the function).
#' @param effects2 An optional factor representing the second treatment effect for two-way designs. Default is NULL.
#' @param effects2.name A character string for the name of the second effect (currently unused in the function).
#' @param k.bonf The Bonferroni correction factor. Must be provided.
#' @param alpha The significance level for the confidence intervals.
#'
#' @return A list containing the following elements:
#'   \item{tau1}{A dataframe of estimated effects for the first treatment.}
#'   \item{tau2}{A dataframe of estimated effects for the second treatment (if applicable).}
#'   \item{general.mean}{A vector of overall means for each response variable.}
#'   \item{p}{The number of response variables (1 for ANOVA, >1 for MANOVA).}
#'   \item{ci.diff.tau1}{A dataframe of confidence intervals for differences between effects of the first treatment.}
#'   \item{ci.diff.tau2}{A dataframe of confidence intervals for differences between effects of the second treatment (if applicable).}
#'   \item{ci.effect.tau1}{A dataframe of confidence intervals for effects of the first treatment.}
#'   \item{ci.effect.tau2}{A dataframe of confidence intervals for effects of the second treatment (if applicable).}
#'
#' @details
#' This function can be used for both MANOVA and ANOVA:
#' - For ANOVA, the input data should have a single response variable (p = 1).
#' - For MANOVA, the input data should have multiple response variables (p > 1).
#' 
#' The function automatically adapts to the number of response variables in the input data.
#' 
#' For one-way designs, use only the `effects1` parameter. For two-way designs, 
#' use both `effects1` and `effects2`.
#'
#' @export
#' @export manova.bonferroni
#'
#' @examples
#' \dontrun{
#' # For one-way ANOVA
#' anova_results <- manova.bonferroni(fit.anova, my_data, effects1 = factor(treatment),
#'                                    k.bonf = 3, alpha = 0.05)
#' 
#' # For one-way MANOVA
#' manova_results <- manova.bonferroni(fit.manova, my_data, effects1 = factor(treatment),
#'                                     k.bonf = 4, alpha = 0.05)
#' 
#' # For two-way ANOVA
#' two_way_anova_results <- manova.bonferroni(fit.anova, my_data, 
#'                                            effects1 = factor(treatment1),
#'                                            effects2 = factor(treatment2),
#'                                            k.bonf = 6, alpha = 0.05)
#' 
#' # For two-way MANOVA
#' two_way_manova_results <- manova.bonferroni(fit.manova, my_data, 
#'                                             effects1 = factor(treatment1),
#'                                             effects2 = factor(treatment2),
#'                                             k.bonf = 8, alpha = 0.05)
#'                                             }
#' 
manova.bonferroni = function(fit.manova, data, effects1, effects2 = NULL, k.bonf, alpha){
  
  if (!is.data.frame(data)) {
    stop("'data' must be a dataframe")
  }
  
  if (!is.factor(effects1)) {
    stop("'effects1' must be a factor")
  }
  
  if (!is.null(effects2) && !is.factor(effects2)) {
    stop("'effects2' must be a factor or NULL")
  }
  
  p = ncol(data)
  
  treatments1 = levels(effects1)
  treatments2 = levels(effects2)
  
  tau1 = matrix(0, nrow = p, ncol = length(treatments1))
  tau1 = data.frame(tau1)
  colnames(tau1) = treatments1
  
  tau2 =  matrix(0, nrow = p, ncol = length(treatments2))
  tau2 = data.frame(tau2)
  colnames(tau2) = treatments2
  
  general.mean = MiniR::colMeans(data)
  
  
  rownames = NULL
  
  for(treat in treatments1){
    rownames = names(tau(data, treat, effects1))
    
    tau1[, treat] = tau(data, treat, effects1)
  }
  
  
  row.names(tau1) = rownames
  rownames = NULL
  # tau1 now contains the effect tau of effect1 
  for(treat in treatments2){
    rownames = names(tau(data, treat, effects2))
    
    tau2[, treat] = tau(data, treat, effects2)
  }
  
  row.names(tau2) = rownames
  
  # now we compute the CI for the differences of effects1
  g = length(treatments1)
  ci.diff = matrix(0, nrow = p*g*(g-1)/2, ncol = 2)
  ci.diff = data.frame(ci.diff)
  colnames(ci.diff) = c("inf", "sup")
  row = 1
  rownames = character(p*g*(g-1)/2)
  for(i in 1:(length(treatments1)-1)){
    for(j in (i+1):length(treatments1))
      for(pos in 1:p){
        # here we compute the CI for the differences between tau1[i], tau1[j] in pos k
        n1 = length(which(effects1 == treatments1[i]))
        n2 = length(which(effects1 == treatments1[j]))
        
        ci.diff[row, ] = CI.diff(tau1[, i], tau1[, j], pos, n1, n2, k.bonf, alpha, fit.manova = fit.manova, p)
        rownames[row] = sprintf("%s - %s [%d]", treatments1[i], treatments1[j], pos)
        row = row + 1
      }
  }
  rownames(ci.diff) = rownames
  ci.diff.tau1 = ci.diff
  
  
  
  
  ci.diff.tau2 = NULL
  if(!is.null(effects2)){
    g = length(treatments2)
    ci.diff = matrix(0, nrow = p*g*(g-1)/2, ncol = 2)
    ci.diff = data.frame(ci.diff)
    colnames(ci.diff) = c("inf", "sup")
    row = 1
    rownames = character(p*g*(g-1)/2)
    for(i in 1:(length(treatments2)-1)){
      for(j in (i+1):length(treatments2))
        for(pos in 1:p){
          # here we compute the CI for the differences between tau1[i], tau1[j] in pos k
          n1 = length(which(effects2 == treatments2[i]))
          n2 = length(which(effects2 == treatments2[j]))
          
          ci.diff[row, ] = CI.diff(tau2[, i], tau2[, j], pos, n1, n2, k.bonf, alpha, fit.manova = fit.manova, p)
          rownames[row] = sprintf("%s - %s [%d]", treatments2[i], treatments2[j], pos)
          row = row + 1
        }
    }
    rownames(ci.diff) = rownames
    ci.diff.tau2 = ci.diff
  }
  
  
  
  
  
  g = length(treatments1)
  ci.effect.tau1 = matrix(0, nrow = g*p, ncol = 2)
  ci.effect.tau1 = data.frame(ci.effect.tau1)
  colnames(ci.effect.tau1) = c("inf", "sup")
  row = 1
  rownames = character(p*g)
  for(i in 1:length(treatments1)){
    for(pos in 1:p){
      n1 = length(which(effects1 == treatments1[i]))
      ci.effect.tau1[row, ] = CI.effect(tau1[,i], pos = pos, n1 = n1, x.mean = general.mean, n.tot = nrow(data),
                                        k.bonf = k.bonf, alpha = alpha, fit.manova = fit.manova, p)
      rownames[row] = sprintf("%s [%d]", treatments1[i], pos)
      row = row + 1                             
      
    }
    
  }
  rownames(ci.effect.tau1) = rownames
  
  ci.effect.tau2 = NULL
  
  if(!is.null(effects2)){
    g = length(treatments2)
    ci.effect.tau2 = matrix(0, nrow = g*p, ncol = 2)
    ci.effect.tau2 = data.frame(ci.effect.tau2)
    colnames(ci.effect.tau2) = c("inf", "sup")
    row = 1
    rownames = character(p*g)
    for(i in 1:length(treatments2)){
      for(pos in 1:p){
        n1 = length(which(effects2 == treatments2[i]))
        ci.effect.tau2[row, ] = CI.effect(tau2[,i], pos = pos, n1 = n1, x.mean = general.mean, n.tot = nrow(data),
                                          k.bonf = k.bonf, alpha = alpha, fit.manova = fit.manova, p)
        rownames[row] = sprintf("%s [%d]", treatments2[i], pos)
        row = row + 1                             
        
      }
      
    }
    rownames(ci.effect.tau2) = rownames
  }
  
  
  
  
  
  result = list(tau1 = tau1, tau2 = tau2, general.mean = general.mean, p = p,
                ci.diff.tau1 = ci.diff.tau1, ci.diff.tau2 = ci.diff.tau2,
                ci.effect.tau1 = ci.effect.tau1, ci.effect.tau2 = ci.effect.tau2)
  return(result)
}



