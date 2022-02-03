source("functions-loading.R")

files <- c("Adult", "Clinical", "Developmental")

# Read special CSV and write it back as normal CSV
str_c("data", files, ".csv") %>%
  map(read.csv.2header) %>%
  setNames(files) -> data

data %>%
  names %>%
  walk(~write_csv(data[[.]], str_c("data", ., ".csv")))

files <- c("Adult", "Clinical", "Developmental", "Genetic")

# Load all files
str_c("data", files, ".csv") %>%
  map(~read.csv(., colClasses = "character", stringsAsFactors = F, check.names = F)) %>%
  setNames(files) -> data

# Load duplicates and test ids
# The duplicates and tests file is a CSV file with three columns:
# - type : either "test" for test account, or "dup" for duplicated accounts"
# - id1 : id for a test account, or old id for duplicated account
# - id2 : new id for duplicated account
read_csv("duplicates_and_tests.csv", col_types = "ccc") %>%
  split(.$type) -> dups_tests

## Remove exact duplicated rows
#data %>%
  #map(distinct) -> data

# Remove test patients
data %>%
  map(~ .x[!.x$Patient.ID %in% dups_tests$test$id1,]) -> data

# Remove duplicated patients
# Comparing all the inputs from all files for all the duplicated pairs, the rows corresponding to id2 contain more data, so this is probably the id in use.
# Let's merge IDs to a common one and let the parsing script delete rows without survey session ids, and choose naturally the most recent data.
# ids are still important to create the link with the genetic data, so we need to keep the same as in the genetic files.
# However, at the moment there is none!!! (except for 2923-3504, labelled as such. So for now I'll use id2 as it is the most recent
# TODO: check the ids in the genetic results file for one of the duplicates

# Rename 2923-3504 to 3504 in the genetic results
data$Genetic %>%
  mutate(Patient.ID = ifelse(Patient.ID == "2923-3504", "3504", Patient.ID)) -> data$Genetic

# Create a list-column dataframes of new names and the levels to collapse
dups_tests$dup %>%
  transmute(names = id2,
            levels = map2(id1, id2, c)) -> dups

# Apply a function that remaps the duplicated entries in every data file
data %>%
  map(function(x)
      {
        # Add the first argument to the list for the parametrized call
        tibble(names = ".f",
               levels = list(x$Patient.ID)) %>%
        bind_rows(dups) -> arguments

        # Do call fct_collapse to map the IDs
        do.call(fct_collapse,
                arguments$levels %>%
                  setNames(arguments$names)) -> newid

        x$Patient.ID <- newid

        x
      }) -> data

# Write back to files
files %>%
  walk(~write_csv(data[[.]], str_c("data", ., ".csv")))
