####################
####################
#' @title convertGoldData
#' @description converts unsigned gold copy data to signed directed adjacency matrix using complete and partial knockdown
#' expression values.  
#' @param goldData .tsv file of all gene interactions but void of signed element
#' @param partialKnockdown dataframe of partial knockdowns of gene expression, wt row is the baseline for comparison
#' @param fullKnockdown dataframe of complete knockdowns of gene expression, wt row is the baseline for comparison
#' @return signed adjacency matrix indicated true connections (activation + and inhibition -)
convertGoldData <- function(goldData, partialKnockdown, fullKnockdown){
  edge_list <-  goldData |> 
    mutate(V1_numeric = as.numeric(gsub("G", "", V1))) |> 
    arrange(V1_numeric) |> 
    select(-V1_numeric)
  
  gene_nums <- sort(as.numeric(gsub("G", "", unique(c(edge_list$V1, edge_list$V2)))))
  genes <- paste0("G", gene_nums)
  
  p <- length(genes)
  signed_adj_matrix <- matrix(0, p, p,
                       dimnames = list(genes, genes))
  
  for (i in 1:nrow(edge_list)) {
    from_gene <- edge_list$V1[i]
    to_gene <- edge_list$V2[i]
    edge <- edge_list$V3[i]
    
    if (edge == 1) {
      baseline_expr     <- fullKnockdown[fullKnockdown$strain == "wt", to_gene]
      partial_removal <- partialKnockdown[partialKnockdown$strain == paste0(from_gene, "(+/-)"), to_gene]
      full_removal    <- fullKnockdown[fullKnockdown$strain == paste0(from_gene, "(-/-)"), to_gene]
     
      delta_partial <- partial_removal - baseline_expr
      delta_full    <- full_removal - baseline_expr 
      
      if (delta_partial < 0 && delta_full < 0) {
        signed_adj_matrix[from_gene, to_gene] <- 1  # activation
      } else if (delta_partial > 0 && delta_full > 0) {
        signed_adj_matrix[from_gene, to_gene] <- -1 # inhibition
      } else {
        if (delta_full < 0) {
          signed_adj_matrix[from_gene, to_gene] <- 1  # activation
        } else if (delta_full > 0) {
          signed_adj_matrix[from_gene, to_gene] <- -1  # inhibition
        } else {
          signed_adj_matrix[from_gene, to_gene] <- 0  # no interaction
        }
      }
    }
  }
  
  return(signed_adj_matrix)
  
  
  
  # compare null mutes and hetData with each connection in the goldData
  # let null mutes be the tie breaker if contradiction exists
  # apply sign to correct connection
  
}
