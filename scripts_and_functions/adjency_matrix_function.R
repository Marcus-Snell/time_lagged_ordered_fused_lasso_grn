library(readr)

# This creates the adjency matrix but...
# I'm not confident it is correct as it semisupervised model will not accept it

# Load datasets
null_mutants <- read_xls("./data/dream2/Network1/InSilico1-null-mutants.xls")
heterozygous <- read_xls("./data/dream2/Network1/InSilico1-heterozygous.xls")

# Extract gene names
genes <- colnames(null_mutants)[-1]  # Assuming first column is 'Time' or WT

# Initialize adjacency matrix (p x p)
p <- length(genes)
adj_matrix <- matrix(0, nrow=p, ncol=p, dimnames=list(genes, genes))

# Populate adjacency matrix based on null-mutant data
for (i in 1:p) {
  for (j in 1:p) {
    if (null_mutants[i, j] < null_mutants[1, j]) {  
      adj_matrix[i, j] <- -1  # Inhibition
    } else if (null_mutants[i, j] > null_mutants[1, j]) {
      adj_matrix[i, j] <- 1  # Activation
    }
  }
}

# Apply heterozygous, **but only if null mutant was inconclusive**
for (i in 1:p) {
  for (j in 1:p) {
    if (adj_matrix[i, j] == 0) {  # Only override if not already set
      if (heterozygous[i, j] < heterozygous[1, j]) {
        adj_matrix[i, j] <- 1  # Activation (heterozygous knockdown)
      } else if (heterozygous[i, j] > heterozygous[1, j]) {
        adj_matrix[i, j] <- -1  # Inhibition
      }
    }
  }
}
