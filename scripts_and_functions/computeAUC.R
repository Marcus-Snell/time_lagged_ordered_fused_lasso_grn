library(pROC)

#' @title computeAUC
#' @description computes the AUC usingn the ground truth data
#' @param adjacencyMatrix matrix of 1s and 0s: output from timeLagLassoNetworkReconstruction()
#' @param goldData vector of true network values, must be preprocessed into a vector where each genes interaction is listed in order
#' @return AUC score and plot of the ROC
computeAUC <- function(adjacencyMatrix = NULL, goldData=NULL) {
  p <- dim(adjacencyMatrix)[1]
  
  predicted_connections <- adjacencyMatrix[!diag(TRUE, p)]
  true_connections <- goldData
  
  roc_curve <- roc(true_connections, predicted_connections)
  
  return(roc_curve)
}
