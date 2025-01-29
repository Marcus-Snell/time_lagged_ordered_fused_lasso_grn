
#data preprocessing functions

#' @title .rescaleData
#' @description .rescaleData
#' @keywords internal
#' @param datasets list of datasets (samples x variables)
#' @return list of datasets with variables scaled by the overall variable means and standard deviations across all input datasets
.rescaleData <- function(datasets){
	# given a list of datasets, combine
	# and rescale each row to mean 0, sd 1
	# and separate
	genes <- colnames(datasets[[1]]) # returns col names, (gene names)
	if(is.list(datasets)){
		datasets <- lapply(datasets, function(ii){ii[,genes]}) # Ensures cols for every dataset in the list contain only contain genes prevoiusly identified
		allData <- do.call(rbind, datasets) #assuming genes x time, combines indvidual datasets into one
		samples <- sapply(datasets, nrow) # vector of the number of rows in each dataset in the list (51, 51)
		samples <- rep(seq_along(samples), samples) # create a seq vector with each num repeated the amount of rows present in a dataframe item
		allData <- scale(allData) # each column values are scaled individually by subtracting the mean of the column and dividing by the SD.
		lapply(unique(samples), function(ii){ # returns the previously combined dataset to the original seperate datasets but with scaled values across all datasets.
			allData[samples==ii, ]
		})
	} else{
		scale(datasets) # scales dataset if there is only one passed (not a list)
	}

#' @title .rescaleDataSeparate
#' @description .rescaleDataSeparate
#' @keywords internal
#' @param datasets list of datasets (samples x variables)
#' @return list of datasets, each separately scaled
.rescaleDataSeparate <- function(datasets){
	# given a list of datasets, treat the datasets separately
	# and rescale each row to mean 0, sd 1
	lapply(datasets, function(ii){ # scales list of datasets as previously stated but with individual mean/sd for each dataset
		scale(ii)
	})
}
