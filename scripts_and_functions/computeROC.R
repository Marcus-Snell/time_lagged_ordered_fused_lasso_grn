library(pROC)

#' @title computeAUC
#' @description computes the AUC usingn the ground truth data
#' @param coeffMatrix matrix of 1s and 0s: output from timeLagLassoNetworkReconstruction()
#' @param goldData vector of true network values, must be preprocessed into a vector where each genes interaction is listed in order
#' @return AUC score and plot of the ROC
computeROC <- function(goldData = NULL, coeffMatrixByLag = NULL) {
  p <- dim(coeffMatrixByLag[[1]])[1]
  
  combined_coeff_matrix <- Reduce('+', coeffMatrixByLag)
  removed_diagonal <- combined_coeff_matrix[!diag(TRUE, p)]
  predicted_connections <- as.vector(removed_diagonal)
  true_connections <- goldData[!diag(TRUE, p)]
  
  roc_curve <- roc(true_connections, predicted_connections)
  
  return(roc_curve)
}
