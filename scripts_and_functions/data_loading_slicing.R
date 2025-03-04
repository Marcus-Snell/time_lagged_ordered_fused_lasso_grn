library(readxl)
library(readr)
library(dplyr)

####################
####################
#' @title loadTsData
#' @description loads DREAM trajectories data (time series) from a .xls or .tsv file
#' @param file_extension accepts a string value either ".xls" or ".tsv"
#' @param file_path accepts a string value to desired directory
#' @return list of time series
loadTsData <- function(file_extension = "", file_path = "") {
  if (file_extension == ".xls") {
    ts_gene_data <- read_xls(file_path)
  
    ts_gene_data <- ts_gene_data |>
      mutate(group = cumsum(is.na(Time))) |>  # Create a group ID based on NA separators
      filter(!is.na(Time))  # Remove NA rows
  
    # Split the dataset into a list of time-series matrices
    ts_gene_data_list <- ts_gene_data |>
      group_by(group) |>
      group_split() |>
      lapply(function(x) x %>% select(-Time, -group))  # Remove 'Time' & 'group' columns
    
  } else {
    ts_gene_data <- read_tsv(file_path)
    
    num_blocks <- length(unique(ts_gene_data$Time))  # Unique time points per experiment
    num_experiments <- nrow(ts_gene_data) / num_blocks  # Total datasets
    
    # Split the dataset into a list
    ts_gene_data_list <- split(ts_gene_data, rep(1:num_experiments, each = num_blocks))
    
    # Remove "Time" column from each matrix
    ts_gene_data_list <- lapply(ts_gene_data_list, function(x) x %>% select(-Time))
  }
  
  return(ts_gene_data_list)
}
 
####################
####################
#' @title loadExprData
#' @description loads DREAM expression data from a .xls or .tsv file
#' @param file_extension accepts a string value either ".xls" or ".tsv"
#' @param file_path accepts a strings value to desired directory
#' @return 
loadExprData <- function(file_extension = "", file_path = "./data/dream2/Network1") {
  if (file_extension == ".xls"){
    files <- list.files(path = file_path, pattern = "*.xls", full.names = TRUE)
    dataset_list <- lapply(files, read_xls)
      
    dataset_list <- lapply(dataset_list, function(df) {
      df <- df[, sapply(df, is.numeric), drop = FALSE]  # Keep only numeric columns, remove "Strain" column
      return(df)
    })
  }
  
  return(dataset_list)
}

####################
####################
#' @title matrix_slice_function
#' @description Loop over data a create datsets with limited time points, start point moves the right(down) each pass
#' @param data_list list of datasets (samples x variables)
#' @param window_size number of time points to be included in each slice
#' @return list of list of datasets with each sublist only containing window_size rows
sliceTsData <- function(data_list, window_size) {
  start <- 1
  stop <- window_size
  reduced_matrix_list = list()
  
  while (nrow(data_list[[1]]) - start >= window_size) {
    temp_list <- list()
    for (df in data_list) {
      small_df <- df[start:stop, ]
      temp_list <- append(temp_list, list(small_df))
    }
    
    reduced_matrix_list <-  append(reduced_matrix_list, list(temp_list))
    start <- start + 1
    stop <- stop + 1
  }
  return(reduced_matrix_list)
}