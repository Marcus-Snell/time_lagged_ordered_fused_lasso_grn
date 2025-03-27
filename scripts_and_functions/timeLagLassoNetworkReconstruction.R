
#' @title .timeLagLassoLaggedData
#' @description .timeLagLassoLaggedData
#' @keywords internal
#' @param xData lagged expression matrix (n x (p * maxLag))
#' @param yData output expression or change-in-expression vector of length n
#' @param maxLag maximum predictor lag
#' @param lambda vector of penalization parameters, containing one or p elements
#' @param intercept if \code{TRUE}, include a model intercept
#' @param beta_pos optional vector containing positive parts of the model coefficients (\code{NULL} or vector of length p * maxLag)
#' @param beta_neg optional vector containing negative parts of the model coefficients (\code{NULL} or vector of length p * maxLag)
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param epsilon convergence error tolerance
#' @return a list of ordered lasso coefficients corresponding to the positive part of the weakly ordered solution (\code{bp}), the negative part of the weakly ordered solution (\code{bn}), the weakly ordered solution (\code{beta}), the weakly ordered intercept (\code{b0}), the strongly ordered intercept (\code{b0.ordered}), and the strongly ordered solution (\code{beta.ordered}). 
.timeLagLassoLaggedData <- function(xData, yData, maxLag, lambda, intercept=TRUE,
                                    beta_pos=NULL, beta_neg=NULL, method="Solve.QP", strongly.ordered=FALSE, 
                                    maxiter=500, inneriter=100, iter.gg=100, epsilon=1e-6){
  
  # a modified version of timeLagLasso that can be used with multiple time series
  # and only returns variables that are (possibly) relevant for network reconstruction
  # assumes data has been normalized/lagged/combined already; x and y have been obtained using functions from dataLagging.R
  
  if (is.null(beta_pos)){
    beta_pos <- rep(0, ncol(xData))
  }
  if (is.null(beta_neg)){
    beta_neg <- rep(0, ncol(xData))
  }
  
  ################
  
  # if one value of lambda is given, use it for all genes and lags;0
  # otherwise, for each gene, use the same penalization across all lags
  if(length(lambda)==1){
    lambda <- rep(lambda, ncol(xData))
  } else{
    lambda <- rep(lambda, each=maxLag)
  }
  
  #solve for coefficients using TL ordered lasso
  est <- .timeLagLassoEstOrdered(x=xData, y=yData, lambda=lambda, maxlag=maxLag, intercept=intercept,
                                 beta_pos=beta_pos, beta_neg=beta_neg,
                                 stdeviation_inverse_scaled=1, standardize=FALSE, strongly.ordered=strongly.ordered,
                                 method=method, maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg, epsilon=epsilon)
  
  list(bp=est$beta_pos, bn=est$beta_neg, beta=est$beta, b0=est$b0,
       b0.ordered=est$b0.ordered, beta.ordered=est$beta.ordered)
}

####################
####################

#' @title .timeLagLassoNetworkLaggedData
#' @description .timeLagLassoNetworkLaggedData
#' @keywords internal
#' @param xData lagged expression matrix (n x (p * maxLag))
#' @param yData output expression or change-in-expression matrix (n x p)
#' @param maxLag maximum predictor lag
#' @param lambda a scalar or p x p matrix of penalization parameters. If a scalar, all coefficients are subject to the same penalization. Otherwise,\code{lambda[j,i]} is the penalization on variable \code{j} in the model for gene \code{i} 
#' @param self if \code{TRUE}, include loops in time-lagged regression models
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered  if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param cores number of parallel cores
#' @import parallel
#' @return a list of coefficient matrices, with each matrix corresponding to a lag and ordered by increasing lag
.timeLagLassoNetworkLaggedData <- function(xData, yData, maxLag, lambda,
                                           self=TRUE, method='Solve.QP', strongly.ordered=TRUE,
                                           maxiter=500, inneriter=100, iter.gg=100,
                                           cores=1){
  
  # assume that x and y have been obtained using functions from dataLagging.R
  
  #if one value of lambda is given, use it for all genes
  if(length(lambda)==1){
    lambda <- matrix(lambda, ncol(yData), ncol(yData))
  }
  
  if(!self){
    # no loops in the regression models
    coeffsByGene <- parallel::mclapply(seq(ncol(yData)), mc.cores=cores, FUN=function(ii){
      # iterate through the genes, learn a model for each
      xInd <- ((ii-1) * maxLag + 1) : (ii * maxLag) #indices corresponding to loops
      .timeLagLassoLaggedData(xData[,-xInd], yData[,ii], maxLag, lambda[-ii, ii],
                              beta_pos=NULL, beta_neg=NULL, method=method, strongly.ordered=strongly.ordered,
                              maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg)
    })
    
    coeffLags <- rep(seq(maxLag), ncol(yData)-1) #lags of each column in y
    
    # create coefficient matrices for each lag and return
    lapply(seq(maxLag), function(ii){
      subsetIndex <- coeffLags == ii
      coeffMatrix <- matrix(0, ncol(yData), ncol(yData),
                            dimnames=list(colnames(yData), colnames(yData)))
      
      #fill in coefficient matrices
      for(jj in seq_along(coeffsByGene)){
        if(strongly.ordered){
          coeffMatrix[-jj,jj] <- coeffsByGene[[jj]]$beta.ordered[subsetIndex]
        } else{
          coeffMatrix[-jj,jj] <- coeffsByGene[[jj]]$beta[subsetIndex]
        }
      }
      
      coeffMatrix
    })
  } else{
    #include loops in the regression models
    coeffsByGene <- parallel::mclapply(seq(ncol(yData)), mc.cores=cores, FUN=function(ii){
      .timeLagLassoLaggedData(xData, yData[,ii], maxLag, lambda[, ii], 
                              beta_pos=NULL, beta_neg=NULL, method=method, strongly.ordered=strongly.ordered,
                              maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg,)
    })
    
    coeffLags <- rep(seq(maxLag), ncol(yData))
    
    # create coefficient matrices for each lag and return
    lapply(seq(maxLag), function(ii){
      subsetIndex <- coeffLags == ii
      coeffMatrix <- matrix(0, ncol(yData), ncol(yData),
                            dimnames=list(colnames(yData), colnames(yData)))
      
      #fill in coefficient matrices
      for(jj in seq_along(coeffsByGene)){
        if(strongly.ordered){
          coeffMatrix[,jj] <- coeffsByGene[[jj]]$beta.ordered[subsetIndex]
        } else{
          coeffMatrix[,jj] <- coeffsByGene[[jj]]$beta[subsetIndex]
        }
      }

      coeffMatrix
    })
  }
}

####################
####################

#' @title .convertCoefficientsToAdjMatrix
#' @description .convertCoefficientsToAdjMatrix
#' @keywords internal
#' @param coeffMatricesByLag list of coefficient matrices, with each matrix corresponding to a lag and ordered by increasing lag
#' @param maxLag maximum lag to use for network prediction, less than or equal to the regression model lag
#' @param epsilon tolerance or threshold for edge prediction
#' @return a predicted adjacency matrix
.convertCoefficientsToAdjMatrix <- function(coeffMatricesByLag, exprDataList, maxLag=1, epsilon=NULL){
  maxLag <- min(maxLag, length(coeffMatricesByLag))
  
  if (is.null(epsilon)) {
    epsilon <- .findOptimalEpsilon(exprDataList)
  }
  
  #threshold absolute value of the coefficients, creates true/false values (binary)
  adjByLag <- lapply(coeffMatricesByLag[1:maxLag], function(ii){
    abs(ii) > epsilon
  })
  
  #aggregate lags for each gene-pair
  adjPredicted <- Reduce('+', adjByLag) > 0
  diag(adjPredicted) <- 0
  
 
  adjPredicted
}

####################
####################

#' @title .findOptimalEpsilon
#' @description .findOptimalEpsilon
#' @keywords internal
#' @param exprDataList list of expression datasets (timepoints x genes) for p genes
#' @return optimal epsilon for coefficient selection
.findOptimalEpsilon <- function(exprDataList, quantileThreshold=0.75) {
  rescaled_data <- .rescaleDataSeparate(exprDataList)  # scale data
  lagged_data <- .transformListMultiLag(rescaled_data) # lag data
  
  lagged_expr_matrix <- lagged_data$xData # select combined lagged predictors
  
  pcor_matrix <- pcor(lagged_expr_matrix)$estimate # find partial correlation matrix
 
  abs_pcor_values <- abs(pcor_matrix[upper.tri(pcor_matrix)])  # Ignore diagonal
  epsilon_optimal <- quantile(abs_pcor_values, quantileThreshold)
  
  return(epsilon_optimal)

}

####################
####################

#combining everything together:
#' @title timeLaggedOrderedLassoNetwork
#' @description Predicts a gene regulatory network from expression data based on the time-lagged ordered lasso.
#' @param exprDataList list of expression datasets (timepoints x genes) for p genes
#' @param output expression output type ('expr' or 'change expr.')
#' @param maxLag maximum predictor lag
#' @param lambda a scalar or p x p matrix of penalization parameters. If a scalar, all coefficients are subject to the same penalization. Otherwise,\code{lambda[j,i]} is the penalization on variable \code{j} in the model for gene \code{i} 
#' @param self if \code{TRUE}, include loops in time-lagged regression models
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param rescale if \code{TRUE}, rescale input exprssion data
#' @param rescaleSeparately if \code{TRUE}, rescale each dataset separately
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param cores number of parallel cores
#' @return a predicted adjacency matrix
#' @export
timeLaggedOrderedLassoNetwork <- function(exprDataList,
                                          output='expr.', maxLag=2, lambda=1, self=TRUE,
                                          method='Solve.QP', strongly.ordered=TRUE,
                                          rescale=TRUE, rescaleSeparately=FALSE, 
                                          maxiter=500, inneriter=100, iter.gg=100,
                                          cores=1){
  
  p <- ncol(exprDataList[[1]])
  
  #checking inputs:
  
  #validate number of columns in each dataset is the same
  if(any(sapply(exprDataList, ncol) != p)){
    stop('exprDataList: must have the same number of genes in each dataset')
  } else if(length(unique(lapply(exprDataList, function(ii){sort(colnames(ii))}))) > 1){
    stop('exprDataList: must have the same genes in each dataset')
  }
  
  # only two types of output accepted, the original gene expression, or the change in expression from the lagged data.
  if(output != 'expr.' && output != 'change expr.'){
    warning("output: must be 'expr.' or 'change expr.'; using 'expr.'", call. = FALSE)
    output <- 'expr.'
  }
  
  # check underlying lasso optimization method: generalized gradient or quadratic programming?
  if(method != 'Solve.QP' && method != 'GG'){
    warning("method: must be 'Solve.QP' or 'GG'; using 'Solve.QP")
    method <- 'Solve.QP'
  }
  
  # coefficient penalty check, single value applys same penalty regardless of time, vector applies different penalty for different lags/coefficients
  if(!is.numeric(lambda) && !is.integer(lambda)){
    stop('lambda: must be of class numeric or integer')
  } else if(length(lambda) != 1 && !all(dim(lambda) == c(p, p))){
    stop('lambda: must be a scalar or p x p matrix')
  }
  
  # max lag check, must be an integer
  if(!is.numeric(maxLag) || maxLag < 1){
    stop('maxLag: must be a scalar greater than or equal to 1')
  } else if(floor(maxLag) < maxLag){
    warning('maxLag: rounding down to an integer')
    maxLag <- floor(maxLag)
  }
  
  # lasso type check, strongly ordered means we think the effects of the coefficients should degrade as time progresses
  if(!is.logical(strongly.ordered)){
    warning('strongly.ordered: must be a logical; using FALSE')
    strongly.ordered <- FALSE
  }
  
  # check if rescale should be done on means and sd of all datasets included
  if(!is.logical(rescale)){
    warning('rescale: must be a logical; using TRUE')
    rescale <- TRUE
    
  }
  
  # check if rescale should be done on means and sd of datasets individually
  if(!is.logical(rescaleSeparately)){
    warning('rescaleSeparately: must be a logical; using FALSE')
    rescaleSeparately <- FALSE
  }
  
  # check the max number of iterations is for convergence is numeric, prevents lasso from running indefinitely
  if(!is.numeric(maxiter) && !is.integer(maxiter)){
    warning('maxiter: must be a scalar; using 500')
    maxiter <- 500
  }
  
  # same as above except this is limiting the number of iteration for each substep
  if(!is.numeric(inneriter) && !is.integer(inneriter)){
    warning('inneriter: must be a scalar; using 100')
    inneriter <- 100
  }
  
  # same as above except controls iterations for when method is generalized gradient (gg)
  if(!is.numeric(iter.gg) && !is.integer(iter.gg)){
    warning('iter.gg: must be a scalar; using 100')
    iter.gg <- 100
  }
  ##########################
  ##########################
  
  #preprocess expr. data, uses dataLagging.R to scale data appropriately
  if(rescale){
    if(rescaleSeparately){
      rescaledList <- .rescaleDataSeparate(exprDataList)
    } else{
      rescaledList <- .rescaleData(exprDataList)
    }
  }
  
  #lag expr. data, uses dataLagging.R to create lagged expression data and the output expression data at time t + 1
  if(output == 'expr.'){
    transformedList <- .transformListMultiLag(rescaledList, maxLag)
  } else if(output == 'change expr.'){
    transformedList <- .transformListMultiLagChange(rescaledList, maxLag)
  }
  
  xData <- transformedList$xData # lagged expression data
  yData <- transformedList$yData # output expression data at time t + 1
  rm(rescaledList, transformedList) # memory cleanup once xData and yData have been created
  
  # creates a p x p matrix to be used as penalization parameter if lambda is a scaler
  if(length(lambda)==1){
    lambda <- matrix(lambda, p, p)
  }
  
  ##########
  
  #compute coefficients
  coefficientsByLag <- .timeLagLassoNetworkLaggedData(xData, yData, maxLag, lambda,
                                                      self=self, method=method, strongly.ordered=strongly.ordered,
                                                      maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg,
                                                      cores=cores)
  
  ##########
  
  #compute adj. matrix
  adjMatrix <- .convertCoefficientsToAdjMatrix(coefficientsByLag, exprDataList, maxLag)
  attr(adjMatrix, "coefficientsByLag") <- coefficientsByLag
  return(adjMatrix)
}



#combining everything together:
#' @title timeLaggedOrderedLassoSemiSupervisedNetwork
#' @description Predicts a posterior gene regulatory network from a prior network and expression data based on the time-lagged ordered lasso.
#' @param exprDataList list of expression datasets (timepoints x genes) for p genes
#' @param adjMatrix prior network adjacency matrix
#' @param output expression output type ('expr' or 'change expr.')
#' @param maxLag maximum predictor lag
#' @param lambdaEdge a scalar edge penalization parameter.
#' @param lambdaNonEdge a scalar non-edge penalization parameter.
#' @param self if \code{TRUE}, include loops in time-lagged regression models
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param rescale if \code{TRUE}, rescale input exprssion data
#' @param rescaleSeparately if \code{TRUE}, rescale each dataset separately
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param cores number of parallel cores
#' @return a predicted adjacency matrix
#' @export
timeLaggedOrderedLassoSemiSupervisedNetwork <- function(exprDataList,
                                                        adjMatrix, 
                                                        output='expr.', maxLag=2, lambdaEdge=1, lambdaNonEdge=1, 
                                                        self=TRUE,
                                                        method='Solve.QP', strongly.ordered=FALSE,
                                                        rescale=TRUE, rescaleSeparately=FALSE, 
                                                        maxiter=500, inneriter=100, iter.gg=100,
                                                        cores=1){
  
  p <- ncol(exprDataList[[1]])
  genes <- colnames(exprDataList[[1]])
  
  #checking inputs:
  
  if(!is.numeric(lambdaEdge) && !is.integer(lambdaEdge)){
    stop('lambda: must be of class numeric or integer')
  } else if(length(lambdaEdge) != 1){
    stop('lambda: must be a scalar')
  }
  
  if(!is.numeric(lambdaNonEdge) && !is.integer(lambdaNonEdge)){
    stop('lambda: must be of class numeric or integer')
  } else if(length(lambdaEdge) != 1){
    stop('lambda: must be a scalar')
  }
  
  if(!all(dim(adjMatrix) == c(p, p))){
    stop('adjMatrix: must be a p x p adjacency matrix')
  }
  adjMatrix <- adjMatrix[genes, genes]
  
  lambda <- adjMatrix * lambdaEdge + (adjMatrix == 0) * lambdaNonEdge
  
  timeLaggedOrderedLassoNetwork(exprDataList,
                                output=output, maxLag=maxLag, lambda=lambda, self=self,
                                method=method, strongly.ordered=strongly.ordered,
                                rescale=rescale, rescaleSeparately=rescaleSeparately, 
                                maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg,
                                cores=cores)
}