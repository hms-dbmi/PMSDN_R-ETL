source("functions-process.R")
library(tidyverse)

export_date = as.Date("2020-12-30")

# Create dir for output, create empty mapping file and ontology object
dir.create("output_transmart", recursive = T)
dir.create("output_human", recursive = T)
cat("Filename\tCategory Code\tColumn Number\tData Label\n", file = "output_transmart/mapping.txt")
ontology <- character(0)

# ==== Demographics ====
processDemographics() %>%
  write_csv("output_human/demographics.csv")

# ==== Phenotypic information ====
processFile("Clinical") %>%
  write_csv("output_human/clinical.csv")
processFile("Adult") %>%
  write_csv("output_human/adult.csv")
processFile("Developmental") %>%
  write_csv("output_human/developmental.csv")

# ==== Genetic information ====
processGenetics() %>%
  write_csv("output_human/genetic.csv")
