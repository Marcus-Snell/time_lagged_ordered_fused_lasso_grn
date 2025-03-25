library(quadprog)
library(ppcor)
library(reshape2)
library(patchwork)

source("./scripts_and_functions/dataPreprocessing.R")
source("./scripts_and_functions/dataLagging.R")
source("./scripts_and_functions/orderedLassoFunctions.R")
source("./scripts_and_functions/timeLagLassoNetworkReconstruction.R")
source("./scripts_and_functions/data_loading_slicing.R")
source("./scripts_and_functions/lambda_selection.R")
source("./scripts_and_functions/plot_functions.R")
source("./scripts_and_functions/computeROC.R")

# Load Data ---------------------------------------------------------------


ts_10_genes <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network2/Time_Series/InSilicoSize10-Ecoli1-trajectories.tsv")
ts_50_genes <- loadTsData(file_extension = ".xls", file_path = "./data/dream2/Network1/Time_series/InSilico1-trajectories.xls")

# True network interactions
ts_10_gold <- read.table("./data/dream3/InSilicoSize10/Network1/gold_copy/DREAM3GoldStandard_InSilicoSize10_Ecoli1.txt")

ts_10_gold <-  ts_10_gold |> 
  mutate(V1_numeric = as.numeric(gsub("G", "", V1))) |> 
  arrange(V1_numeric) |> 
  select(-V1_numeric)

true_vals <- 
  ts_10_gold |> 
  select(V3) |> 
  as.vector()

true_vals <- c(true_vals[[1]])
  

# Time Lagged Lasso Network Reconstruction --------------------------------


# Find optimal lamdas matrix
opt_lambdas_10_genes <- findOptimalLambdas(ts_10_genes, combine_data = FALSE)
opt_lambdas_50_genes <- findOptimalLambdas(ts_50_genes, num_folds = 10)

# Creates predicted(inferred) GRN based on time series data of different gene expression levels
time_lag_ord_lasso_net_10_genes <- timeLaggedOrderedLassoNetwork(ts_10_genes, lambda = opt_lambdas_10_genes,  maxLag = 2)
time_lag_ord_lasso_net_50_genes <- timeLaggedOrderedLassoNetwork(ts_50_genes, lambda = opt_lambdas_50_genes)


# Setup data for gene network evolution visualization. Creates a list of adjacency matrices that can be plotted
adj_matrix_list <- list()
sliced_ts_data_list <- sliceTsData(ts_10_genes, 10) # slice data

# Create adj matrices for each sliced dataset list and add to adj_matrix_list
opt_lambdas_list <- list()

# Find optimal lambdas for time series slices
for (list in sliced_ts_data_list) {
  temp_lambda <- findOptimalLambdas(list)
  opt_lambdas_list <- append(opt_lambdas_list, list(temp_lambda))
}

# Run time lagged ordered lasso on all slices with found optimal lambdas
for (i in seq_along(sliced_ts_data_list)) {
    temp_adj <- timeLaggedOrderedLassoNetwork(sliced_ts_data_list[[i]], lambda = opt_lambdas_list[[i]])
    adj_matrix_list <- append(adj_matrix_list, list(temp_adj))
}



# Evaluation Metrics ------------------------------------------------------


# Plot ROC, find AUC and plot 
coeffs_by_lag <- attr(time_lag_ord_lasso_net_10_genes, "coefficientsByLag")
roc_curve <- computeROC(true_vals, coeffs_by_lag)
auc(roc_curve)
plot(roc_curve)

# Find ROC, AUC for all slices of time series
roc_auc_list <- list()
for (i in seq_along(adj_matrix_list)) {
  coeffs_by_lag <- attr(adj_matrix_list[[i]], "coefficientsByLag")
  roc_curve <- computeROC(true_vals, coeffs_by_lag)
  auc_score <- auc(roc_curve)
  roc_auc_list <- append(roc_auc_list, list(auc_score, roc_curve))
}

# plot ROCS for small time series
roc_list <- roc_auc_list[seq(2, length(roc_auc_list), by = 2)]
auc_only_list <- roc_auc_list[seq(1, length(roc_auc_list), by = 2)]

roc_plot_list <- list()
ts_num <- 1
for (i in seq_along(roc_list)) {
  roc <- roc_list[[i]]
  auc_val <- auc_only_list[[i]] 
  df <- coords(roc, "all", transpose = FALSE)
  
  x_pos <- max(1 - df$specificity, na.rm = TRUE) * 0.80
  y_pos <- min(df$sensitivity, na.rm = TRUE) + 0.15
  
  roc_plot <- 
    ggplot(df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "blue", linewidth = 1) +
    annotate("text", x = x_pos, y = y_pos, label = "AUC: ") +
    annotate("text", x = x_pos + 0.1, y = y_pos, label = round(auc_val, 3), color = "red") + 
    labs(x = "False Positive Rate",
         y = "True Positive Rate",
         caption = paste0("GRN Evolution ", ts_num))
  
  roc_plot_list <- append(roc_plot_list, list(roc_plot))
  ts_num <- ts_num + 1
}

small_ts_roc <- wrap_plots(roc_plot_list, ncol = 3) +
                plot_annotation(title = "GRN Evolution ROC Curves", 
                                subtitle = "DREAM3 InSilico - 10 Gene Network",
                theme = theme(plot.title = element_text(size = 18),
                              plot.subtitle = element_text(size = 12),
                              axis.title.x = element_text(size = 2, face = "bold"),
                              axis.title.y = element_text(size = 2, face = "bold")))

small_ts_roc

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
