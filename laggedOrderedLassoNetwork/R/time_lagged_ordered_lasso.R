source("dataPreprocessing.R")
source("dataLagging.R")
source("orderedLassoFunctions.R")
source("timeLagLassoNetworkReconstruction.R")

library(readxl)
library(quadprog)
library(igraph)
library(ggraph)

# Processing Data Files ---------------------------------------------------


# Import gene expression data
files <- list.files(path = "C:/Users/Marcu/Documents/Data_Science/Spring 2025/Capstone II/laggedOrderedLassoNetwork/Data/dream2/Network1", pattern = "*.xls", full.names = TRUE)
dataset_list <- lapply(files, read_xls)

dataset_list <- lapply(dataset_list, function(df) {
  df <- df[, sapply(df, is.numeric), drop = FALSE]  # Keep only numeric columns, remove "Strain" column, is this right?
  return(df)
})

# Import trajectory data
adj_matrix <- read_xls("C:\\Users\\Marcu\\Documents\\Data_Science\\Spring 2025\\Capstone II\\laggedOrderedLassoNetwork\\Data\\dream2\\InSilico1-trajectories.xls")


# Execute Time Lag Lasso Network Reconstruction ---------------------------


# the lines below work but where does the trajectories.xls data come in to play?
time_lag_ord_lasso_net <- timeLaggedOrderedLassoNetwork(dataset_list)
time_lag_ord_lasso_net_test <- timeLaggedOrderedLassoNetwork(dataset_list, output = "change expr.")


# outputs matrix of 1s and 0s.  this represents connectivity between genes
timeLaggedOrderedLassoSemiSupervisedNetwork(dataset_list, time_lag_ord_lasso_net)


# Create Plots ------------------------------------------------------------


# create directed graph plot
testg <- graph_from_adjacency_matrix(time_lag_ord_lasso_net)
plot(testg)

#another package to experiment with for directed graph
graph <- as_tbl_graph(time_lag_ord_lasso_net) |> 
  mutate(Popularity = centrality_degree(mode = 'in'))

ggraph(time_lag_ord_lasso_net, layout = 'fr') +
  geom_edge_link() +
  geom_node_point(size = 9, color = "dodgerblue3") +
  theme_graph(background = 'white') + 
  geom_node_text(label = colnames(time_lag_ord_lasso_net), color = "white", size = 3)
