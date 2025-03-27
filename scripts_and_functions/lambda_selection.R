library(glmnet)
#library(orderedLasso) # use instead of glmnet, will have to downgrade R version

#' @title lambda_selection
#' @description Loop over data to find optimal lambda for each gene value, uses glmnet cross validation method
#' @param expr_data_list list of time series expression data
#' @param min minimum value for lambda
#' @param max maximum value for lambda
#' @param by amount of increase for each lambda as it advance to the max value
#' @param num_folds number of folds for cross validation
#' @return matrix of optimal lambdas for use in timeLaggedOrderedLassoNetwork
findOptimalLambdas <- function(expr_data_list, min=0, max=10, by=0.1, num_folds=5, combine_data = TRUE){
  possible_lambdas <- seq(min, max, by) # create possible lambdas
  
  p <- ncol(expr_data_list[[1]])  # find col dimensions for initial matrix and col/row dimensions for final matrix
  x <- length(expr_data_list) # find row dimension for initial separate matrices
  
  #scale data
  rescaled_data_list <- .rescaleData(expr_data_list)
  
  if (combine_data) {
    combined_data <- do.call(rbind, rescaled_data_list)
    
    optimal_lambdas_vector <- numeric(p)  # store optimal lambda for each gene
    
    for (j in 1:p) {
      response <- combined_data[, j]                    # one gene as response
      predictors <- combined_data[, -j, drop = FALSE]   # all other genes as predictors
      
      lasso_cv <- cv.glmnet(predictors, response, alpha = 1, nfolds = num_folds)
      optimal_lambdas_vector[j] <- lasso_cv$lambda.min
    }
    
    # Create a matrix where each optimal lambda fills its entire column
    optimal_lambdas_matrix <- matrix(rep(optimal_lambdas_vector, each = p), nrow = p, ncol = p)
    
  } else {
  
  #find optimal lambda using glmnet(), possibly orderedLasso
  temp_lambda_matrix <- matrix(0, x, p)
  for (i in seq_along(rescaled_data_list)) {
    data_list <- rescaled_data_list[[i]]
    for (j in 1:p) {
      response <- data_list[, j] # select one column
      predictors <- data_list[, -j, drop = FALSE] # select all other columns
    
      # find optimal lambda for a response then add it to the temp matrix
      lasso_cross_validation <- cv.glmnet(predictors, response, alpha = 1, nfolds = num_folds) # conduct cross validation
      best_lambda <- lasso_cross_validation$lambda.min # find optimal lambda for gene
      temp_lambda_matrix[i, j] <- best_lambda # add optimal to temp matrix
      }
  }
  
  # average columns to find optimal lambda per gene
  optimal_lambdas_vector <- colMeans(temp_lambda_matrix)
  optimal_lambdas_matrix <- matrix(rep(optimal_lambdas_vector, each = p), p, p)
  
    }
  
  return(optimal_lambdas_matrix)
}