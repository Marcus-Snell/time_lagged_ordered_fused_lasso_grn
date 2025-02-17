source("./scripts_and_functions/dataPreprocessing.R")
source("./scripts_and_functions/dataLagging.R")
source("./scripts_and_functions/orderedLassoFunctions.R")
source("./scripts_and_functions/timeLagLassoNetworkReconstruction.R")

library(readxl)
library(quadprog)
library(dplyr)
library(readr)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(cowplot)

# Processing Data Files ---------------------------------------------------


# Import gene expression data (this is the ground truth expression data).  This is used for the adj matrix
# that is passed to the semisupervised model (null-mutants and heterozygous data).
files <- list.files(path = "./data/dream2/Network1", pattern = "*.xls", full.names = TRUE)
dataset_list <- lapply(files, read_xls)

dataset_list <- lapply(dataset_list, function(df) {
  df <- df[, sapply(df, is.numeric), drop = FALSE]  # Keep only numeric columns, remove "Strain" column, is this right?
  return(df)
})

# Import trajectory data (time series data), this is what we will pass to timeLaggedOrderedLassoNetwork()
# Two seperate gene networks to work with here
ts_data_50_genes <- read_xls("./data/dream2/Network1/Time_series/InSilico1-trajectories.xls")

ts_data_50_genes <- ts_data_50_genes %>%
  mutate(group = cumsum(is.na(Time))) %>%  # Create a group ID based on NA separators
  filter(!is.na(Time))  # Remove NA rows

# Split the dataset into a list of time-series matrices
ts_data_list_50_genes <- ts_data_50_genes %>%
  group_by(group) %>%
  group_split() %>%
  lapply(function(x) x %>% select(-Time, -group))  # Remove 'Time' & 'group' columns


ts_data_10_genes <- read_tsv("./data/dream3/InSilicoSize10/Network1/Time_Series/InSilicoSize10-Ecoli1-trajectories.tsv")

num_blocks <- length(unique(ts_data_10_genes$Time))  # Unique time points per experiment
num_experiments <- nrow(ts_data_10_genes) / num_blocks  # Total datasets

# Split the dataset into a list
ts_data_list_10_genes <- split(ts_data_10_genes, rep(1:num_experiments, each=num_blocks))

# Remove "Time" column from each matrix
ts_data_list_10_genes <- lapply(ts_data_list_10_genes, function(x) x %>% select(-Time))

# Execute Time Lag Lasso Network Reconstruction ---------------------------


# Creates predicted(inferred) GRN based on time series data of different gene expression levels
time_lag_ord_lasso_net_50_genes <- timeLaggedOrderedLassoNetwork(ts_data_list_50_genes)
time_lag_ord_lasso_net_10_genes <- timeLaggedOrderedLassoNetwork(ts_data_list_10_genes)

# Utilizes same time series data as well as previously know GRN info to output the predicted GRN
#timeLaggedOrderedLassoSemiSupervisedNetwork(ts_data_list, adj_matrix)


# Slice Gene Expression Data and Create Smaller Dataframes ----------------------


#' @title .matrix_slice_function
#' @description Loop over data a create datsets with limited time points, start point moves the right(down) each pass
#' @param data_list list of datasets (samples x variables)
#' @param window_size number of time points to be included in each slice
#' @return list of list of datasets with each sublist only containing window_size rows
matrix_slice_function <- function(data_list, window_size) {
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

# creates list of adjacency matrices that can be plotted
adj_matrix_list <- list()
sliced_data_list <- matrix_slice_function(ts_data_list_10_genes, 10) # slice data

# create adj matrices for each slice dataset list and add to adj_matrix_list
for (list in sliced_data_list) {
  temp_adj <- timeLaggedOrderedLassoNetwork(list)
  adj_matrix_list <- append(adj_matrix_list, list(temp_adj))
  
}


# Create Plots ------------------------------------------------------------


# create directed graph plot
grn_plot_function <- function(grn) {
  graph <- as_tbl_graph(grn, directed = TRUE) %>%
    activate(edges) %>%
    mutate(
      interaction_type = ifelse(weight > 0, "activator", "inhibitor"),  # Identify activation/inhibition
      edge_color = ifelse(weight > 0, "dodgerblue", "pink2"),   # Green for activation, red for inhibition
      edge_linetype = ifelse(weight > 0, "solid", "solid") # Solid for activation, dashed for inhibition
    )

  # these lines aren't doing anything just yet.  Need to figure out how to add two geom_edge_link()
  activator_arrow <- arrow(length = unit(3, 'mm'), type = "closed")  # Standard arrow
  inhibitor_arrow <- arrow(length = unit(3, 'mm'), type = "open", ends = "last")

  # plot graph
  ggraph(graph, layout = 'linear', circular = TRUE) +
    geom_edge_link(aes(color = interaction_type),
                   arrow = arrow(length = unit(3, 'mm')),
                   start_cap = circle(5, 'mm'),
                   end_cap = circle(5, 'mm')
                   ) +
    geom_node_point(size = 9, color = "yellow3") +
    theme_graph(background = 'white') +
    geom_node_text(aes(label = name), color = "midnightblue", size = 3) + 
    scale_edge_color_manual(
      values = c("activator" = "dodgerblue", "inhibitor" = "pink2"),  # Assign colors properly
      name = "Interaction Type"
    ) +
    theme_graph(background = "grey13") + 
    theme(legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"))
  }

# plot single network
full_10_gene_plot <- grn_plot_function(time_lag_ord_lasso_net_10_genes) + 
  labs(title = "Gene Expression and Repression") + 
  theme(plot.title = element_text(color = "white"))

ggsave(filename = "10_gene_grn_plot.png", plot = full_10_gene_plot, width = 6.5, height = 5)

# plot sliced data and facet, needs work still the scale isn't correct
small_grn_plots <- lapply(adj_matrix_list, grn_plot_function)

# save small GRN plots individually
index <- 1
for (i in seq_along(small_grn_plots)) {
  modified_plot <- small_grn_plots[[i]] + 
    labs(title = paste0("GRN Evolution ", index)) +  # Dynamically generate title
    theme(plot.title = element_text(color = "white"))
  
  ggsave(filename = paste0("small_plot_", index, ".png"), plot = modified_plot, width = 6.5, height = 5, dpi = 300)
  index <- index + 1
}


# group all plots together on one chart
final_plot <- plot_grid(plotlist = small_grn_plots, ncol = 4, align = "hv", 
                        rel_widths = rep(1, length(small_grn_plots)), 
                        rel_heights = rep(1, length(small_grn_plots)))
  

print(final_plot)

ggsave("GRN_plots.png", final_plot, width = 25, height = 12, dpi = 300)
