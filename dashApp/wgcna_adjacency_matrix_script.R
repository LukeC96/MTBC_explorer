calculate_adjacency_matrix <- function(){
  # Check if WGCNA is installed, if not, install it
  if (!requireNamespace("WGCNA", quietly = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("WGCNA")
  }

  library(WGCNA)

  require("knitr")
  knitr::opts_chunk$set(echo = TRUE)
  # Set the root/top-level directory to the project directory so all paths
  # that follow will be relative to that directory
  opts_knit$set(root.dir = "../")
  require("here")

  # expression counts for genes from each sample
  analysis <- readRDS(here("R_data/datExpr.RData"))

  # create correlation adjacency matrix each gene to each gene (or use distance)
  adjMat <- adjacency(analysis, # nolint: object_name_linter.
    selectCols = NULL,
    type = "unsigned",
    power = 14,
    corFnc = "bicor", corOptions = "maxPOutliers=0.05",
    weights = NULL
  )

  return(adjMat)
}