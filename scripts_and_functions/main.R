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
source("./scripts_and_functions/convertGoldData.R")


# Load Data ---------------------------------------------------------------


# gene expression time series
ecoli1_ts <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network1/Time_Series/InSilicoSize10-Ecoli1-trajectories.tsv")
ecoli2_ts <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network2/Time_Series/InSilicoSize10-Ecoli2-trajectories.tsv")

yeast1_ts <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network3/Time_Series/InSilicoSize10-Yeast1-trajectories.tsv")
yeast2_ts <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network4/Time_Series/InSilicoSize10-Yeast2-trajectories.tsv")
yeast3_ts <- loadTsData(file_extension = ".tsv", file_path = "./data/dream3/InSilicoSize10/Network5/Time_Series/InSilicoSize10-Yeast3-trajectories.tsv")


# True network interactions
ecoli1_gold <- read.table("./data/dream3/InSilicoSize10/Network1/gold_copy/DREAM3GoldStandard_InSilicoSize10_Ecoli1.txt")
ecoli1_partial_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network1/InSilicoSize10-Ecoli1-heterozygous.tsv")
ecoli1_full_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network1/InSilicoSize10-Ecoli1-null-mutants.tsv")

ecoli1_adj <- convertGoldData(ecoli1_gold, ecoli1_partial_knockdown, ecoli1_full_knockdown)
ecoli1_true_vals <- as.vector(ecoli1_adj)

ecoli2_gold <- read.table("./data/dream3/InSilicoSize10/Network2/gold_copy/DREAM3GoldStandard_InSilicoSize10_Ecoli2.txt")
ecoli2_partial_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network2/InSilicoSize10-Ecoli2-heterozygous.tsv")
ecoli2_full_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network2/InSilicoSize10-Ecoli2-null-mutants.tsv")

ecoli2_adj <- convertGoldData(ecoli2_gold, ecoli2_partial_knockdown, ecoli2_full_knockdown)
ecoli2_true_vals <- as.vector(ecoli2_adj)


yeast1_gold <- read.table("./data/dream3/InSilicoSize10/Network3/gold_copy/DREAM3GoldStandard_InSilicoSize10_Yeast1.txt")
yeast1_partial_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network3/InSilicoSize10-Yeast1-heterozygous.tsv")
yeast1_full_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network3/InSilicoSize10-Yeast1-null-mutants.tsv")

yeast1_adj <- convertGoldData(yeast1_gold, yeast1_partial_knockdown, yeast1_full_knockdown)
yeast1_true_vals <- as.vector(yeast1_adj)

yeast2_gold <- read.table("./data/dream3/InSilicoSize10/Network4/gold_copy/DREAM3GoldStandard_InSilicoSize10_Yeast2.txt")
yeast2_partial_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network4/InSilicoSize10-Yeast2-heterozygous.tsv")
yeast2_full_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network4/InSilicoSize10-Yeast2-null-mutants.tsv")

yeast2_adj <- convertGoldData(yeast2_gold, yeast2_partial_knockdown, yeast2_full_knockdown)
yeast2_true_vals <- as.vector(yeast2_adj)

yeast3_gold <- read.table("./data/dream3/InSilicoSize10/Network5/gold_copy/DREAM3GoldStandard_InSilicoSize10_Yeast3.txt")
yeast3_partial_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network5/InSilicoSize10-Yeast3-heterozygous.tsv")
yeast3_full_knockdown <- read_tsv("./data/dream3/InSilicoSize10/Network5/InSilicoSize10-Yeast3-null-mutants.tsv")

yeast3_adj <- convertGoldData(yeast3_gold, yeast3_partial_knockdown, yeast3_full_knockdown)
yeast3_true_vals <- as.vector(yeast3_adj)
  

# Time Lagged Lasso Network Reconstruction --------------------------------


# Find optimal lamdas matrix
ecoli1_opt_lambdas <- findOptimalLambdas(ecoli1_ts, combine_data = FALSE)
ecoli2_opt_lambdas <- findOptimalLambdas(ecoli1_ts, combine_data = FALSE)

yeast1_opt_lambdas <- findOptimalLambdas(yeast1_ts, combine_data = FALSE)
yeast2_opt_lambdas <- findOptimalLambdas(yeast2_ts, combine_data = FALSE)
yeast3_opt_lambdas <- findOptimalLambdas(yeast2_ts, combine_data = FALSE)


# Creates predicted(inferred) GRN based on time series data of different gene expression levels
ecoli1_grn <- timeLaggedOrderedLassoNetwork(ecoli1_ts, lambda = ecoli1_opt_lambdas, maxLag = 1)
ecoli2_grn <- timeLaggedOrderedLassoNetwork(ecoli2_ts, lambda = ecoli2_opt_lambdas, maxLag = 1)

yeast1_grn <- timeLaggedOrderedLassoNetwork(yeast1_ts, lambda = yeast1_opt_lambdas, maxLag = 1)
yeast2_grn <- timeLaggedOrderedLassoNetwork(yeast2_ts, lambda = yeast2_opt_lambdas, maxLag = 1)
yeast3_grn <- timeLaggedOrderedLassoNetwork(yeast3_ts, lambda = yeast3_opt_lambdas, maxLag = 1)


# Run gene network evolution. Creates a list of adjacency matrices that can be plotted
grnEvolution <- function(timeSeries) {
  adj_matrix_list <- list()
  opt_lambdas_list <- list()
  sliced_ts_data_list <- sliceTsData(timeSeries, 10) # slice data

  # Find optimal lambdas for time series slices
  for (ts_slice in sliced_ts_data_list) {
    temp_lambda <- findOptimalLambdas(ts_slice)
    opt_lambdas_list <- append(opt_lambdas_list, list(temp_lambda))
  }

  # Run time lagged ordered lasso on all slices with found optimal lambdas
  for (i in seq_along(sliced_ts_data_list)) {
      temp_adj <- timeLaggedOrderedLassoNetwork(sliced_ts_data_list[[i]], lambda = opt_lambdas_list[[i]])
      adj_matrix_list <- append(adj_matrix_list, list(temp_adj))
  }
  return(adj_matrix_list)
}

ecoli1_adj_list <- grnEvolution(ecoli1_ts)
ecoli2_adj_list <- grnEvolution(ecoli2_ts)

yeast1_adj_list <- grnEvolution(yeast1_ts)
yeast2_adj_list <- grnEvolution(yeast2_ts)
yeast3_adj_list <- grnEvolution(yeast3_ts)


# Evaluation Metrics ------------------------------------------------------


# Plot ROC, find AUC and plot 
coeffs_by_lag <- attr(yeast3_grn, "coefficientsByLag")
roc_curve <- computeROC(yeast3_true_vals, coeffs_by_lag)
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

small_ts_roc_plot <- wrap_plots(roc_plot_list, ncol = 3) +
                plot_annotation(title = "GRN Evolution ROC Curves", 
                                subtitle = "DREAM3 InSilico - 10 Gene Network",
                theme = theme(plot.title = element_text(size = 18),
                              plot.subtitle = element_text(size = 12),
                              axis.title.x = element_text(size = 2, face = "bold"),
                              axis.title.y = element_text(size = 2, face = "bold")))

ggsave(path = "./visualizations", filename = "small_ts_auc_scores.png", plot = small_ts_roc_plot, width = 14.5, height = 10)

# Visualize Gene Networks -------------------------------------------------


# plot single networks
full_10_gene_plot <- plotGrnDirectedGraph(time_lag_ord_lasso_net_10_genes) + 
  labs(title = "Gene Expression and Repression - 10 Gene Network") + 
  theme(plot.title = element_text(color = "white", size = 13))

ggsave(path = "./visualizations", filename = "10_gene_grn_plot.png", plot = full_10_gene_plot, width = 6.5, height = 5)

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
