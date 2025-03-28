library(ggplot2)
library(ggraph)
library(tidygraph)
library(cowplot)
library(ggdark)

####################
####################
#' @title plotGrnDirectedGraph
#' @description plots directed graph using timeLaggedOrderedLassoNetwork() output 
#' @param grn adjacency matrix
#' @return ggplot chart of GRN
plotGrnDirectedGraph <- function(grn) {
  graph <- as_tbl_graph(grn, directed = TRUE) |>
  activate(edges) |>
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

####################
####################
#' @title plotaucLambdaLagComparison
#' @description plots directed graph using timeLaggedOrderedLassoNetwork() output 
#' @param grnResultsDF a dataframe of different AUC scores run on different lambda, and lag values
#' @param networkName a string value of the name of the network plotted (ex. "Ecoli1")
#' @return ggplot chart of individual lambdas and their respective scores at each lag value
plotaucLambdaLagComparison <- function(grnResultsDF, networkName){
  grnResultsDF |>
  ggplot(aes(x = max_lag, y = auc, color = as.factor(lambda), group = lambda)) + 
    geom_line() + 
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(color = "Lambda",
         title = "Model Performance Comparison by Varying Max Lag and Lambda",
         subtitle = paste0("Network: ", networkName)) +
    dark_theme_light() +
    theme(plot.title = element_text(size = 22),
          plot.subtitle = element_text(size = 18),
          legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
}