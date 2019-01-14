source("functions-loading.R")
library(readxl)

 files <- c("Clinical", "Developmental", "Adult")

# Transform HTML files to special CSV
str_c("data", files) %>%
  walk(html2csv)

# Transform genetic results from XLSX to CSV
read_xlsx("dataGenetic.xlsx") %>%
  rename(Patient.ID = ID) %>%
  filter(!Patient.ID %>% is.na) %>%
  write_csv("dataGenetic.csv")
