library(quadprog)

source("./scripts_and_functions/dataPreprocessing.R")
source("./scripts_and_functions/dataLagging.R")
source("./scripts_and_functions/orderedLassoFunctions.R")
source("./scripts_and_functions/timeLagLassoNetworkReconstruction.R")
source("./scripts_and_functions/data_loading_slicing.R")
source("./scripts_and_functions/lambda_selection.R")
source("./scripts_and_functions/plot_functions.R")


# Load Data ---------------------------------------------------------------


ts_10_genes <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network1/Time_Series/InSilicoSize10-Ecoli1-trajectories.tsv")
ts_50_genes <- loadTsData(file_extension = ".xls", file_path = "./data/dream2/Network1/Time_series/InSilico1-trajectories.xls")


# Time Lagged Lasso Network Reconstruction --------------------------------


# Find optimal lamdas matrix
opt_lambdas_10_genes <- cv_lambda_selection(ts_10_genes)
opt_lambdas_50_genes <- cv_lambda_selection(ts_50_genes)

# Creates predicted(inferred) GRN based on time series data of different gene expression levels
time_lag_ord_lasso_net_10_genes <- timeLaggedOrderedLassoNetwork(ts_10_genes, lambda = opt_lambdas_10_genes)
time_lag_ord_lasso_net_50_genes <- timeLaggedOrderedLassoNetwork(ts_50_genes, lambda = opt_lambdas_50_genes)

# Setup data for gene network evolution visualization. Creates a list of adjacency matrices that can be plotted
adj_matrix_list <- list()
sliced_ts_data_list <- sliceTsData(ts_10_genes, 10) # slice data

# create adj matrices for each sliced dataset list and add to adj_matrix_list
#TODO: Write a loop to find optimal lambdas for each sliced dataset
for (list in sliced_ts_data_list) {
  temp_adj <- timeLaggedOrderedLassoNetwork(list)
  adj_matrix_list <- append(adj_matrix_list, list(temp_adj))
  
}


# Visualize Gene Networks -------------------------------------------------


# plot single networks
full_10_gene_plot <- plotGrnDirectedGraph(time_lag_ord_lasso_net_10_genes) + 
  labs(title = "Gene Expression and Repression - 10 Gene Network") + 
  theme(plot.title = element_text(color = "white", size = 13))

ggsave(path = "./visualizations", filename = "10_gene_grn_plot.png", plot = full_10_gene_plot, width = 6.5, height = 5)

full_50_gene_plot <- plotGrnDirectedGraph(time_lag_ord_lasso_net_50_genes) + 
  labs(title = "Gene Expression and Repression - 50 Gene Network") + 
  theme(plot.title = element_text(color = "white", size = 25))

ggsave(path = "./visualizations", filename = "50_gene_grn_plot.png", plot = full_50_gene_plot, width = 11.5, height = 10)

# Plot individual networks from sliced data and save (allows visualization of Gene Network Evolution)
small_grn_plots <- lapply(adj_matrix_list, plotGrnDirectedGraph)

index <- 1
for (i in seq_along(small_grn_plots)) {
  modified_plot <- small_grn_plots[[i]] + 
    labs(title = paste0("GRN Evolution ", index)) +  # Dynamically generate title
    theme(plot.title = element_text(color = "white"))
  
  ggsave(path = "./visualizations", filename = paste0("small_plot_", index, ".png"), plot = modified_plot, width = 6.5, height = 5, dpi = 300)
  index <- index + 1
}

# Group all GRN evolution plots together on one chart
final_plot <- plot_grid(plotlist = small_grn_plots, ncol = 4, align = "hv", 
                        rel_widths = rep(1, length(small_grn_plots)), 
                        rel_heights = rep(1, length(small_grn_plots)))


ggsave(path = "./visualizations", "all_grn_plots.png", final_plot, width = 25, height = 12, dpi = 300)
