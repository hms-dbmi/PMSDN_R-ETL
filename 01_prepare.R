source("functions-loading.R")
library(readxl)

 files <- c("Clinical", "Developmental", "Adult")

# Transform HTML files to special CSV
str_c("data", files) %>%
  walk(html2csv)

# Transform genetic results from XLSX to CSV
#
# A couple things to do first :
#
# Rename the file from Catalina as dataGenetic.xlsx
# (-> or better, get Calatina to name it that way when she sends it to you)
#
# Add any additional case to this file, and check that tho columns match
# (-> or get Catalina to send a complete file, keeping with the same structure (no additional column, no renaming))
#
# Patient ID 10370 has multiple misplaced cells : Array2_Year, Array2_Lab, Array2_Method, Array2_Version, CodingDNAChange, ProteinChange, and Effect. For now this is censored.
# (-> get Catalina to fix this in her file)
# Mark Patient 2679 with an "Incomplete" Test verification, and enquire for the correct data

read_xlsx("dataGenetic.xlsx") %>%
  rename(Patient.ID = ID) %>%
  filter(!Patient.ID %>% is.na) %>%
  write_csv("dataGenetic.csv")
