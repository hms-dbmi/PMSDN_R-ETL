source("functions-mapping.R")
source("functions-reformatting.R")

# Process at the file level
processFile <- function(questionnaire)
{
  # Add the questionnaire level to the ontology
  ontology <<- push(ontology, questionnaire)

  # Read the data and premapping files
  data <- read.csv(paste0("data", questionnaire, ".csv"), colClasses = "character", stringsAsFactors = F, check.names = F)

  cat(paste0("└┬ Questionnaire : ", questionnaire, "\n"))

  data <- data[data$Survey.Session.ID != "", ]

  premap <- read.csv(paste0("premap", questionnaire, ".csv"), stringsAsFactors = F, colClasses = "character")

  premap$ColNum <- as.integer(premap$ColNum)

  data2 <- data["Patient.ID"] %>% distinct()

  # Process each SubFile level (excluding the empty SubFile level->Demographics)
  for (subfile in levels(factor(premap$SubFile, exclude = "")))
    data2 <- merge(data2, processSubfile(questionnaire, subfile, data, premap), by = "Patient.ID")

  ontology <<- pop(ontology)

  data2
}

# Process at the SubFile level
processSubfile <- function(questionnaire, subfile, data, premap)
{
  cat(paste0(" └┬ Subfile : ", subfile, "\n"))
  # Add the SubFile level to the ontology
  ontology <<- push(ontology, subfile)

  # Subset the premapping file with only the current SubFile
  premap <- filter(premap, SubFile == subfile)

  # Create new data frame to contain transformed/curated data
  data2 <- data["Patient.ID"] %>% distinct()

  # Process each Head1 level and merge resulting data
  for (head1 in unique(premap$Head1))
    data2 <- merge(data2, processHead1(head1, data, premap), by = "Patient.ID")

  # Parse resulting var names and write mappings
  addMappings(questionnaire, subfile, ontology, data2)

  # Write the $SubFile.txt
  write.table(data2, file = paste0("output_transmart/", questionnaire, "-", subfile, ".txt"), row.names = F, sep = "\t", quote = F, na = "")

  ontology <<- pop(ontology)

  data2
}

# Process at the Head1 level
processHead1 <- function(head1, data, premap)
{
  cat(paste0("  ├─ Header 1 : ", head1, "\n"))
  # Subset the premapping file with only the current Head1
  premap <- filter(premap, Head1 == head1)

  # Anchor-based filtering of variables from the data file
  data <- anchorFilter(premap, data)

  # Sort by Survey.Time
  data$Survey.Time <- as.numeric(strptime(data$Survey.Time, format = "%Y-%m-%d %H:%M:%S"))
  data <- arrange(data, Patient.ID, Survey.Time)

  # Delete records made less than 24 hours before the next
  data <- filter(data, (lead(Survey.Time) - Survey.Time) > 24 * 3600 | lead(Patient.ID) != Patient.ID | Patient.ID == max(Patient.ID))

 # Reformatting needed. Execute the function given in the premap file for the reformatting
  if (any(premap$Reformat != ""))
  {
    funcname <- levels(factor(premap$Reformat, exclude = ""))
    cat(paste0("  │             + executing function ", funcname, "\n"))
    eval(parse(text = paste0("data <- ", funcname, "(data, premap)")))
  }

  # Manage 'Other Value' columns
  if (any(grepl("_Other.Value$", names(data))))
  {
    cat("  │             + replacing Other.Value\n")
    data <- otherValue(data)
  }

  # Manage "checkbox" items
  if (any(grepl("_Unsure$", names(data))))
  {
    cat("  │             + processing checkboxes\n")
    data <- checkboxes(data)
  }

  # Manage longitudinal data
  if (any(premap$Evo == "1"))
  {
    cat("  │             + processing evolutive data\n")
    data <- evolutive(data)
  }
  else
    data <- historical(data)

  # Clean the content of tabs and linefeeds
  for (var in names(data))
    data[var] <- gsub("[\n\t]", " ", data[[var]])

  data
}

processDemographics <- function()
{
  # Read raw data files
  adult         <- read.csv("dataAdult.csv", colClasses = "character", check.names = F, stringsAsFactors = F, na.strings = "")
  developmental <- read.csv("dataDevelopmental.csv", colClasses = "character", check.names = F, stringsAsFactors = F, na.strings = "")
  clinical      <- read.csv("dataClinical.csv", colClasses = "character", check.names = F, stringsAsFactors = F, na.strings = "")

  # Delete rows with no Survey Session ID
  adult         <- adult[adult$Survey.Session.ID != "", ]
  developmental <- developmental[developmental$Survey.Session.ID != "", ]
  clinical      <- clinical[clinical$Survey.Session.ID != "", ]

  # Extract basic demographic informations (patient ID, SEX, AGE, RACE, COUNTRY)
  adult[                    c("Patient.ID", "Birthdate", "Gender", "Race", "Country")] %>%
    bind_rows(clinical[     c("Patient.ID", "Birthdate", "Gender", "Race", "Country")]) %>%
    bind_rows(developmental[c("Patient.ID", "Birthdate", "Gender", "Race", "Country")]) %>%
    unique() %>%
    arrange(Patient.ID) %>%
    group_by(Patient.ID) %>%
    mutate(Birthdate = as.Date(Birthdate)) %>%
    mutate(Age = as.numeric(export_date - Birthdate) / 365.25) %>%
    mutate(Age_months = Age * 12) -> Demographics

  # Write Demographics.txt
  write.table(Demographics, "output_transmart/Demographics.txt", row.names = F, sep = "\t", quote = F)

  # Write the mappings
  ontology <<- push(ontology, "Demographics")
  addMapping("Demographics.txt", ontology, 1, "SUBJ_ID")
  addMapping("Demographics.txt", ontology, 2, "BIRTHDATE")
  addMapping("Demographics.txt", ontology, 3, "SEX")
  addMapping("Demographics.txt", ontology, 4, "RACE")
  addMapping("Demographics.txt", ontology, 5, "COUNTRY")
  addMapping("Demographics.txt", ontology, 6, "AGE_IN_YEARS")
  addMapping("Demographics.txt", ontology, 7, "AGE")
  ontology <<- pop(ontology)

  Demographics
}

# Prepare names for reshaping and transmart mapping
rename_patient <- function(x)
{
  x %>%
    str_replace(., "_", "\\.") %>% 
    str_replace("(\\w+)_([\\w.]+)_(\\d)", "\\1_\\3_\\2") 
}

processGenetics <- function()
{
  ontology <<- push(ontology, "Genetic")

  # Read the curated genetic data
  genetic <- read.csv("dataGenetic.csv", stringsAsFactors = F, na.strings = "NA")

  # ==== Overall genetic testing status ====
  genetic %>%
    select(Patient.ID,
           GeneticStatus,
           SHANK3Involved,
           PathogenicSHANK3Defect,
           Mosaicism,
           Origin,
           TestVerification) -> geneticStatus

  write.table(geneticStatus, file = paste0("output_transmart/", "Genetics-Status.txt"), row.names = F, sep = "\t", quote = F, na = "")
  addMappings("Genetics", "Overall Status", ontology, geneticStatus)

  # Sequencing results
  genetic %>%
    filter(Result.type == "Sequencing") %>%
    select(Patient.ID,
           # Result.type,
           SequencingEffect,
           SequencingClinicalSignificance) %>%
    set_names(c("Patient.ID", "Sequencing_Effect", "Sequencing_Clinical.Significance")) -> genetic_sequencing

  #Array results
  genetic %>%
    filter(Result.type == "Array") %>%
    select(Patient.ID, starts_with("Array")) %>%
    reshape(direction = "wide",
            idvar = "Patient.ID",
            timevar = "Array.num",
            sep = "_") -> genetic_temp

  genetic_array_name <- rename_patient(colnames(genetic_temp))
  genetic_temp %>%
    set_names(genetic_array_name) -> genetic_array

  genetic_sequencing %>%
    full_join(genetic_array) -> genetics


  write.table(genetics, file = paste0("output_transmart/","Genetics-Results.txt"), row.names = F, sep = "\t", quote = F, na = "")
  addMappings("Genetics","Results",ontology,genetics)
  ontology <<- pop(ontology)

  genetic
}
