source("functions-loading.R")
source("functions-mapping.R")
source("functions-reformatting.R")
source("functions-genes.R")

# Process at the file level
processFile <- function(questionnaire, noOutput = F)
{
  # Add the questionnaire level to the ontology
  if (!noOutput)
    ontology <<- push(ontology, questionnaire)

  # Read the data and premapping files
  print(paste("processing data for ", questionnaire))
  data <- read.csv.2header(paste0("data", questionnaire, ".csv"))
    print("Finished parsing data")

  data <- data[data$Survey.Session.ID != "", ]

  print(paste("processing premap for ", questionnaire))
  premap <- read.csv(paste0("premap", questionnaire, ".csv"), stringsAsFactors = F, colClasses = "character")
      print("Finished premap")
    
  premap$ColNum <- as.integer(premap$ColNum)

  data2 <- data["Patient.ID"] %>% distinct()

  # Process each SubFile level (excluding the empty SubFile level->Demographics)
  for (subfile in levels(factor(premap$SubFile, exclude = "")))
      print(paste("Processing ", subfile))
      data2 <- merge(data2, processSubfile(questionnaire, subfile, data, premap, noOutput = noOutput), by = "Patient.ID")
      print("Finished subfile")

  if (!noOutput)
    ontology <<- pop(ontology)

  if (noOutput)
    return(data2)
}

# Process at the SubFile level
processSubfile <- function(questionnaire, subfile, data, premap, noOutput)
{
  # Add the SubFile level to the ontology
  if (!noOutput)
    ontology <<- push(ontology, subfile)

  # Subset the premapping file with only the current SubFile
  premap <- filter(premap, SubFile == subfile)

  # Create new data frame to contain transformed/curated data
  data2 <- data["Patient.ID"] %>% distinct()

  # Process each Head1 level and merge resulting data
  for (head1 in unique(premap$Head1))
    data2 <- merge(data2, processHead1(head1, data, premap), by = "Patient.ID")

  # Parse resulting var names and write mappings
  if (!noOutput)
  {
    addMappings(questionnaire, subfile, ontology, data2)

    # Write the $SubFile.txt
    write.table(data2, file = paste0("output/", questionnaire, "-", subfile, ".txt"), row.names = F, sep = "\t", quote = F, na = "")

    ontology <<- pop(ontology)
  }

  data2
}

# Process at the Head1 level
processHead1 <- function(head1, data, premap)
{
  # Subset the premapping file with only the current Head1
  premap <- filter(premap, Head1 == head1)

  # Anchor-based filtering of variables from the data file
  data <- anchorFilter(premap, data)

  # Sort by Survey Date
  data$Survey.Date <- as.numeric(strptime(data$Survey.Date, format = "%Y-%m-%d %H:%M:%S"))
  data <- arrange(data, Patient.ID, Survey.Date)

  # Delete records made less than 24 hours before the next
  data <- filter(data, (lead(Survey.Date) - Survey.Date) > 24 * 3600 | lead(Patient.ID) != Patient.ID | Patient.ID == max(Patient.ID))

 # Reformatting needed. Execute the function given in the premap file for the reformatting
  if (any(premap$Reformat != ""))
  {
    funcname <- levels(factor(premap$Reformat, exclude = ""))
    eval(parse(text = paste0("data <- ", funcname, "(data, premap)")))
  }

  # Manage 'Other Value' columns
  if (any(grepl("_Other.Value$", names(data))))
    data <- otherValue(data)

  # Manage "checkbox" items
  if (any(grepl("_Unsure$", names(data))))
    data <- checkboxes(data)

  # Manage longitudinal data
  if (any(premap$Evo == "1"))
    data <- evolutive(data)
  else
    data <- historical(data)
  
  # Clean the content of tabs and linefeeds
  for (var in names(data))
    data[var] <- gsub("[\n\t]", " ", data[[var]])

  data
}

processDemographics <- function(noOutput = F)
{
  # Read raw data files
  adult         <- read.csv.2header("dataAdult.csv")
  developmental <- read.csv.2header("dataDevelopmental.csv")
  clinical      <- read.csv.2header("dataClinical.csv")

  # Delete rows with no Survey Session ID
  adult         <- adult[adult$Survey.Session.ID != "", ]
  developmental <- developmental[developmental$Survey.Session.ID != "", ]
  clinical      <- clinical[clinical$Survey.Session.ID != "", ]

  # Extract basic demographic informations (patient ID, SEX, AGE, RACE, COUNTRY)
  adult[                    c("Patient.ID", "Birthdate", "Gender", "Ancestral.Background", "Country")] %>%
    bind_rows(clinical[     c("Patient.ID", "Birthdate", "Gender", "Ancestral.Background", "Country")]) %>%
    bind_rows(developmental[c("Patient.ID", "Birthdate", "Gender", "Ancestral.Background", "Country")]) %>%
    unique() %>%
    arrange(Patient.ID) %>%
    group_by(Patient.ID) %>%
    mutate(Birthdate = as.Date(Birthdate)) %>%
    mutate(Age = as.numeric(export_date - Birthdate) / 365.25) %>%
    mutate(Age_months = Age * 12) -> Demographics

  if (!noOutput)
  {
    # Write Demographics.txt
    write.table(Demographics, "output/Demographics.txt", row.names = F, sep = "\t", quote = F)

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
  }

  if (noOutput)
    return(Demographics)
}

processGenetics <- function(noOutput = F)
{
  if (!noOutput)
    ontology <<- push(ontology, "Genetic")

  # Read the curated genetic data
  Genetics <- read.csv("dataGenetics.csv", stringsAsFactors = F, na.strings = "")

  # ==== "Misc" genetic data ====
  if (!noOutput)
  {
    ontology <<- push(ontology, "Misc")
      Genetics_misc <- select(Genetics,
                              Patient.ID,
                              Test.Date,
                              Comments,
                              Std.Nomenclature,
                              Test.Method,
                              Gene,
                              Laboratory,
                              Category,
                              Test.Verification,
                              Consultant.Verification,
                              Karyotype.Start,
                              Karyotype.End,
                              Array.Version,
                              Array.Confirmation.Studies,
                              Parental.Results,
                              Parental.Test.Method,
                              Parental.Origin)
      Genetics_misc$Comments <- gsub("[\n\t]", " ", Genetics_misc$Comments)
      Genetics_misc$Std.Nomenclature <- gsub("[\n\t]", " ", Genetics_misc$Std.Nomenclature)
  
      write.table(Genetics_misc, file = paste0("output/", "Genetics-Misc.txt"), row.names = F, sep = "\t", quote = F, na = "")
      addMappings("Genetics", "Misc", ontology, Genetics_misc)
    ontology <<- pop(ontology)
  }
    
  # ==== "Alterations" genetic data ====
  if (!noOutput)
    ontology <<- push(ontology, "Chromosomal alterations")
  
  # Reshape the data frame to prepare it for the coordinates conversion
  varying = grep("\\d$", names(Genetics), value = T)
  Genetics <- reshape(Genetics, direction = "long", varying = varying, sep = ".", idvar = "Patient.ID", timevar = "Test.nb") %>%
    arrange(Patient.ID, Test.nb) %>%
    filter(!is.na(Gain_Loss))
  
  Genetics_ranges   <- processRanges(Genetics)
  
  # Keep only GRCh38/hg38 deletion ranges
  Genetics_ranges <- filter(Genetics_ranges, Genome.Browser.Build == "GRCh38/hg38") %>% select(-Genome.Browser.Build)
  
  if (noOutput)
  {
    return(Genetics_ranges)
  }
  else
  {
    # Reshape the table in the wide format
    Genetics_ranges <- group_by(Genetics_ranges, Patient.ID, Chr_Gene) %>% mutate(N = row_number()) %>% ungroup() %>% ungroup()
    
    data2 <- distinct(Genetics_ranges[c("Patient.ID", "Result.type")])
    for (chr in unique(Genetics_ranges$Chr_Gene))
    {
      ranges <- filter(Genetics_ranges, Chr_Gene == chr) %>% select(-Chr_Gene, -Result.type)
      for (n in unique(ranges$N))
      {
        ranges2 <- filter(ranges, N == n)
        names(ranges2) <- c("Patient.ID", paste0("Chr", chr, "_", n,"_",c("Gain/Loss","Start","End")), "N")
        data2 <- merge(data2, select(ranges2, -N), by = "Patient.ID", all = T)
      }
    }
    
    write.table(data2, file = paste0("output/","Genetics-Ranges.txt"), row.names = F, sep = "\t", quote = F, na = "")
    addMappings("Genetics","Ranges",ontology,data2)
    ontology <<- pop(ontology)
  
    ontology <<- pop(ontology)
  }
}
