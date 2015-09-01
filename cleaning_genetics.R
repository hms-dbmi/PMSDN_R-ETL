source("functions-process.R")
source("functions-mapping.R")

# Read the curated genetic data
Genetics <- read.csv("dataGenetics.csv", stringsAsFactors = F, na.strings = "")

# ==== Genes ====
Genetics_genes    <- processGenes(Genetics)

# ==== Pathways ====
Genetics_pathways <- processPathways(Genetics, Genetics_genes)
