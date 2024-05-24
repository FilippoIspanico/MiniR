#' Plot Confidence Intervals
#'
#' This function plots Bonferroni and/or Simultaneous T2 confidence intervals for a given set of data.
#'
#' @param BonfCI A matrix or vector representing the Bonferroni confidence intervals. If a vector is provided, it will be reshaped into a matrix with 3 columns.
#' @param SimCI A matrix or vector representing the Simultaneous T2 confidence intervals. If a vector is provided, it will be reshaped into a matrix with 3 columns.
#' @param delta.0 A vector representing the hypothesized values under H0.
#' @param legend.position A character string indicating the position of the legend. Default is 'bottomright'.
#'
#' @details
#' The function plots the confidence intervals provided. If both `BonfCI` and `SimCI` are given, both sets of intervals will be plotted. If only one is provided, only that set of intervals will be plotted.
#' 
#' - The Bonferroni intervals are plotted in orange.
#' - The Simultaneous T2 intervals are plotted in blue.
#' - The hypothesized values under H0 are plotted in black.
#'
#' @return A plot showing the confidence intervals.
#'
#' @examples
#' # Example usage with both Bonferroni and Simultaneous T2 intervals
#' BonfCI <- matrix(c(1, 2, 3, 2, 3, 4, 1, 3, 4), ncol = 3, byrow = TRUE)
#' SimCI <- matrix(c(1.5, 2.5, 3.5, 2.5, 3.5, 4.5, 1.5, 3.5, 4.5), ncol = 3, byrow = TRUE)
#' delta.0 <- c(2, 3, 3)
#' plot.intervals(BonfCI, SimCI, delta.0)
#'
#' @export

plot.intervals <- function(BonfCI = NULL, SimCI = NULL, delta.0, legend.position = 'bottomright') {
  
  
  if(is.null(BonfCI) && is.null(SimCI)){
    stop("Both BonfCI and SimCI cannot be NULL")
  }
  
  if(!is.null(BonfCI)){
    
    if (is.null(dim(BonfCI))) {
      BonfCI <- matrix(BonfCI, ncol = 3, byrow = TRUE)
    }
    
    n_intervals = nrow(BonfCI)
    
    matplot(t(matrix(1:n_intervals, n_intervals, 3)), t(BonfCI), type='b', pch='', xlim=c(1-.05,n_intervals+0.5), xlab='',
            ylab='', main='Confidence intervals')
    
    # Plotting the Bonferroni intervals
    segments(matrix(1:n_intervals, n_intervals, 1), BonfCI[,1], matrix(1:n_intervals, n_intervals, 1), BonfCI[,3],
             col='orange', lwd=2)
    points(1:n_intervals, BonfCI[,2], col='orange', pch=16)
    
    # Plotting delta.0 under H0 
    points(1:n_intervals+.02, delta.0, col='black', pch=16)
    
    # Plotting the simultaneous T2
    if(is.null(SimCI)){
      
      legend(legend.position, c('Bonf. IC'), col=c('orange'), lty=1, lwd=2)
      return()
    }
    
    segments(matrix(1:n_intervals+.1,n_intervals,1),SimCI[,1],matrix(1:n_intervals+.1,n_intervals,1),SimCI[,3], col='blue', lwd=2)
    points(1:n_intervals+.1,SimCI[,2], col='blue', pch=16)
    
    legend(legend.position, c('Bonf. IC', 'Sim-T2 IC'), col=c('orange', 'blue'), lty=1, lwd=2)
    
  }
  
  else{
    
    if (is.null(dim(SimCI))) {
      SimCI <- matrix(SimCI, ncol = 3, byrow = TRUE)
    }
    
    n_intervals = nrow(SimCI)
    
    
    matplot(t(matrix(1:n_intervals, n_intervals, 3)), t(SimCI), type='b', pch='', xlim=c(1-.05,n_intervals+0.5), xlab='',
            ylab='', main='Confidence intervals')
    segments(matrix(1:n_intervals,n_intervals,1),SimCI[,1],matrix(1:n_intervals,n_intervals,1),SimCI[,3], col='blue', lwd=2)
    points(1:n_intervals, SimCI[,2], col='blue', pch=16)
    
    legend(legend.position, c('Sim-T2 CI'), col=c('blue'), lty=1, lwd=2)
    
    points(1:n_intervals, delta.0, col='black', pch=16)
    
  }
}




#' Interpret Principal Component Analysis Results
#'
#' This function generates a comprehensive set of plots to aid in the interpretation of principal component analysis (PCA) results. It produces three main plots: the variance of the principal components, the variance of the original variables, and the cumulative contribution of the principal components to the total variance. Additionally, it provides barplots of the loadings for a specified number of principal components.
#'
#' @param pc.result A list object containing the PCA result, typically obtained from the \code{\link{prcomp}} or \code{\link{princomp}} functions. It should include at least \code{$sdev} and \code{$loadings} components.
#' @param data A data frame or matrix containing the original data that was used for the PCA. This data is used to compute the variances of the original variables.
#' @param number_pcs An integer indicating the number of principal components to plot the loadings for. The default is 3.
#'
#' @return This function does not return any value. It generates plots as a side effect.
#'
#' @details The function generates the following plots:
#' \itemize{
#'   \item A scree plot showing the variances of the principal components.
#'   \item A bar plot showing the variances of the original variables.
#'   \item A line plot showing the cumulative contribution of the principal components to the total variance.
#'   \item Bar plots showing the loadings for each of the specified principal components.
#' }
#' 
#' @note The function modifies the graphical parameters using \code{par} to arrange multiple plots on one page. It restores the original parameters after plotting.
#'
#' @examples
#' \dontrun{
#' data(iris)
#' pca_result <- prcomp(iris[, -5], scale. = TRUE)
#' pc.interpretation(pca_result, iris[, -5])
#' }
#'
#' @importFrom graphics layout plot barplot abline box axis par grid
#'
#' @export
#' 

pc.interpretation = function (pc.result, data, number_pcs = 3) {
  
  layout(matrix(c(2, 3, 1, 3), 2, byrow = TRUE))
  plot(pc.result, las = 2, main = 'Principal components')
  
  barplot(sapply(data, var), las = 2, main = 'Original Variables',
          ylab = 'Variances')
  
  plot(cumsum(pc.result$sdev^2) / sum(pc.result$sdev^2), type = 'b', axes = FALSE,
       xlab = 'Number of components', ylab = 'Contribution to the total variance', ylim = c(0, 1))
  
  abline(h = 1, col = 'blue')
  abline(h = 0.8, lty = 2, col = 'blue')
  box()
  axis(2, at = 0:10/10, labels = 0:10/10)
  axis(1, at = 1:ncol(data), labels = 1:ncol(data), las = 2)
  par(mfrow=c(1,1))
  
  
  # Loadings
  loadings.data <- pc.result$loadings
  
  par(mfrow=c(number_pcs,1))
  for(i in 1:number_pcs){
    barplot(loadings.data[,i], ylim = c(-1, 1), main=paste('Loadings PC ',i,sep=''))
    grid()
  }
  par(mfrow=c(1,1))
}




#' Plot Decision Boundaries for LDA and QDA Models
#'
#' This function plots the decision boundaries for Linear Discriminant Analysis (LDA) and Quadratic Discriminant Analysis (QDA) models using ggplot2. It visualizes the separation of different classes in the feature space based on the model's predictions.
#'
#' @param formula A formula specifying the model. The left-hand side of the formula should be the response variable, and the right-hand side should be the predictor variables.
#' @param data A data frame containing the variables in the formula.
#' @param fit.result A fitted model object of class \code{lda} or \code{qda} from the \code{\link[MASS]{lda}} or \code{\link[MASS]{qda}} functions, respectively.
#'
#' @return A ggplot object displaying the decision boundaries and the data points. Points are colored by their true class, and the background is shaded according to the predicted class.
#'
#' @details The function works as follows:
#' \itemize{
#'   \item It extracts the feature variables and response variable from the formula and data.
#'   \item It creates a grid of points spanning the feature space.
#'   \item It predicts the class for each point in the grid using the provided LDA or QDA model.
#'   \item It uses ggplot2 to plot the decision boundaries and data points.
#' }
#' 
#' @note The function currently supports only two feature variables.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' library(ggplot2)
#' library(gridExtra)
#' 
#' data(iris)
#' lda_model <- lda(Species ~ Sepal.Length + Sepal.Width, data = iris)
#' plot_decision_boundaries(Species ~ Sepal.Length + Sepal.Width, iris, lda_model)
#' }
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_tile scale_color_manual scale_fill_manual labs theme_minimal ggtitle
#' @importFrom stats model.frame model.matrix model.response
#'
#' @export
#' 
#' 

plot_decision_boundaries <- function(formula, data, fit.result) {
  # Extract features and labels from the formula and data
  mf <- model.frame(formula, data)
  features <- model.matrix(formula, mf)[, -1]  # Remove intercept column
  labels <- model.response(mf)
  
  feature_names <- colnames(features)
  
  # Create a grid of points spanning the feature space
  x_min <- min(features[, 1]) - 1
  x_max <- max(features[, 1]) + 1
  y_min <- min(features[, 2]) - 1
  y_max <- max(features[, 2]) + 1
  grid <- expand.grid(V1 = seq(x_min, x_max, length.out = 200), 
                      V2 = seq(y_min, y_max, length.out = 200))
  
  # Ensure the grid has the correct column names
  colnames(grid) <- colnames(features)
  
  # Predict the class for each point in the grid
  if ("lda" %in% class(fit.result)) {
    grid$pred <- predict(fit.result, grid)$class
    model_name <- "LDA"
  } else if ("qda" %in% class(fit.result)) {
    grid$pred <- predict(fit.result, grid)$class
    model_name <- "QDA"
  } else {
    stop("The fit.result object must be of class 'lda' or 'qda'")
  }
  
  # Convert the features and labels into a data frame for ggplot
  data_plot <- data.frame(features, label = labels)
  
  # Plot the decision boundary using ggplot2
  p <- ggplot(data_plot, aes_string(x = feature_names[1], y = feature_names[2])) +
    geom_point(aes(color = factor(label)), size = 3) +
    geom_tile(data = grid, aes_string(x = feature_names[1], y = feature_names[2], fill = 'factor(pred)'), alpha = 0.3) +
    scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +
    scale_fill_manual(values = c("red", "blue", "green", "purple", "orange")) +
    labs(color = "Class", fill = "Class", x = feature_names[1], y = feature_names[2]) +
    theme_minimal() +
    ggtitle(paste(model_name, "Decision Boundaries"))
  
  # Return the plot
  return(p)
}



