
#functions for lagging data and constructing expression and change-in-expression output variables

#' @title .transformMultiLag
#' @description .transformMultiLag
#' @keywords internal
#' @param exprData gene expression matrix (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output expression data at time t+1, respectively
.transformMultiLag <- function(exprData, maxLag=1){
	# expression at time t+1
	yData <- exprData[-(1:maxLag),] # create yData, remove the first 1-maxLag rows and all columns
	k <- nrow(yData)  # get number of rows in yData, original number rows minus maxlag rows from beginning

	if(k > 0){
		genes <- colnames(exprData) # get the column names for dataset

		# lag xData
		if(maxLag > 1){
			xData <- lapply(maxLag:1, function(ii){ # creates new columns depending on lag #, if lag == 2, each gene will now have 2 col each
				laggedMatrix <- exprData[(ii):(ii+k-1),] # decrement column lags starting with 1 and moving to maxlag
				colnames(laggedMatrix) <- paste(genes,maxLag-ii+1,sep='-') # rename columns adding "-" and lag num from 1 to max lag
				laggedMatrix
			})
			xData <- do.call(cbind, xData)  # combine columns into 1 datafram (xData)
		} else{
			xData <- exprData[1:k,] # xData will be the same as orginal minus the first row
			colnames(xData) <- paste(genes, 1, sep='-')  # rename columns adding "-" and lag num (only 1 in this instance)
		}
		#reorder so the each gene's lags are blocked together and ordered by increasing lag
		reordered <- paste(rep(genes, each=maxLag), rep(seq(maxLag), ncol(exprData)), sep='-')
		xData <- xData[,reordered]

		lags <- rep(1:maxLag, length(genes)) # ??
		genes <- rep(genes, maxLag) # ??

		rownames(yData) <- NULL  # removes row names for yData
		rownames(xData) <- NULL  # removes row names for xData

		list(xData=xData, yData=yData)  # create list with both data frames, xData and yData
	} else{
		list(xData=NULL, yData=NULL)
	}

}

#' @title .transformListMultiLag
#' @description .transformListMultiLag
#' @keywords internal
#' @param exprDataList list of gene expression matrices (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output expression data at time t+1, respectively
.transformListMultiLag <- function(exprDataList, maxLag=1){
	#apply transformMultiLag to a list of expr. datasets with the same set of genes in the same order
	genes <- colnames(exprDataList[[1]]) # returns col names, (gene names)
	transformed <- lapply(exprDataList, function(curExpr){ # uses transformMultiLag to create xdata and ydata where ydata is the lagged matrix (xdata - maxlag)
		.transformMultiLag(curExpr[, genes], maxLag)
	})

	list(xData=do.call(rbind, lapply(transformed, function(ii){ii$xData})),
		yData=do.call(rbind, lapply(transformed, function(ii){ii$yData})))
}

#' @title .transformMultiLagChange
#' @description .transformMultiLagChange
#' @keywords internal
#' @param exprData gene expression matrix (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data;xData and the output change-in-expression data at time t+1, respectively
.transformMultiLagChange <- function(exprData, maxLag=1){
	# change in expression at time t+1
	yData <- exprData[(maxLag+1):nrow(exprData),] - exprData[maxLag:(nrow(exprData)-1),]
	k <- nrow(yData)

	if(k > 0){
		genes <- colnames(exprData)

		# lag xData
		if(maxLag > 1){
			xData <- lapply(maxLag:1, function(ii){
				laggedMatrix <- exprData[(ii):(ii+k-1),]
				colnames(laggedMatrix) <- paste(genes,maxLag-ii+1,sep='-')
				laggedMatrix
			})
			xData <- do.call(cbind, xData)
		} else{
			xData <- exprData[1:k,]
			colnames(xData) <- paste(genes, 1, sep='-')
		}
		#reorder so the each gene's lags are blocked together and ordered by increasing lag
		reordered <- paste(rep(genes, each=maxLag), rep(seq(maxLag), ncol(exprData)), sep='-')
		xData <- xData[,reordered]

		lags <- rep(1:maxLag, length(genes))
		genes <- rep(genes, maxLag)

		rownames(yData) <- NULL
		rownames(xData) <- NULL

		list(xData=xData, yData=yData)
	} else{
		list(xData=NULL, yData=NULL)
	}

}

#' @title .transformListMultiLagChange
#' @description .transformListMultiLagChange
#' @keywords internal
#' @param exprDataList list of gene expression matrices (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output change-in-expression data at time t+1, respectively
.transformListMultiLagChange <- function(exprDataList, maxLag=1){
	#apply transformMultiLagChange to a list of expr. datasets with the same set of genes in the same order
	genes <- colnames(exprDataList[[1]])
	transformed <- lapply(exprDataList, function(curExpr){
		.transformMultiLagChange(curExpr[, genes], maxLag)
	})

	list(xData=do.call(rbind, lapply(transformed, function(ii){ii$xData})),
		yData=do.call(rbind, lapply(transformed, function(ii){ii$yData})))
}
