source("functions-loading.R")
library(lubridate)
library(magrittr)

# Extract the genetic test results fields from the clinical file, only where there are results
read.csv.2header("dataGenetic.csv") %>%
  select(Patient.ID,
         Test.Date,
         Genetic.Status,
         Genomic.Position.Start.Min,
         Genomic.Position.End.Min,
         Genomic.Position.Start.Max,
         Genomic.Position.End.Max,
         Participant.Origin,
         Comments,
         Std.Nomenclature,
         Test.Method,
         Genome.Browser.Build,
         Gene,
         Lab,
         Category,
         Test.Verification,
         `Results.Verified.-.Consultant`,
         Karyotype.Start,
         Karyotype.End,
         Array.Version,
         Array.Confirmation.Studies,
         `Start.Exon/Intron`:Parental.Origin,
         `Other.Gain/Loss.1`,
         `Other.Arm.Gain/Loss.1`,
         Other.Position.Start.1,
         Other.Position.End.1,
         Origin.1,
         `Other.Gain/Loss.2`,
         `Other.Arm.Gain/Loss.2`,
         Other.Position.Start.2,
         Other.Position.End.2,
         Origin.2,
         `Other.Gain/Loss.3`,
         `Other.Arm.Gain/Loss.3`,
         Other.Position.Start.3,
         Other.Position.End.3,
         Origin.3) %>%
  filter(Genetic.Status != "No Results Received") %>%
# Parse the dates
  mutate(Test.Date = gsub("/", "-",                      Test.Date, perl = T)) %>%
  mutate(Test.Date = gsub("^(\\d)-", "0\\1-",            Test.Date, perl = T)) %>%
  mutate(Test.Date = gsub("-(\\d)-", "-0\\1-",           Test.Date, perl = T)) %>%
  mutate(Test.Date = gsub("-(\\d{2})$", "-20\\1",        Test.Date, perl = T)) %>%
  mutate(Test.Date = gsub("(^\\d{2}-\\d{4}$)", "01-\\1", Test.Date, perl = T)) %>%
  mutate(Test.Date = gsub("^[\\w ]+$", "",               Test.Date, perl = T)) %>%
  mutate(Test.Date = mdy(Test.Date)) %>%
  distinct -> Genetics

# Correct variable names
names(Genetics) <- gsub("/", ".", names(Genetics))

# Correct test dates in the 20th century
Genetics$Test.Date[which(Genetics$Test.Date > as.Date(now()))] <- Genetics$Test.Date[which(Genetics$Test.Date > as.Date(now()))] - years(100)

# Sort by patient and date
Genetics %<>%
  arrange(Patient.ID, Test.Date) %>%
# Create new variables with defaults and reorder the columns
  mutate(Result.type = "coordinates",
         Gain_Loss.1 = "Loss",
         Chr_Gene.1 = "22",
         Start.1 = Genomic.Position.Start.Min,
         End.1 = Genomic.Position.End.Min,
         Origin.1 = Participant.Origin,
         Gain_Loss.2 = Other.Gain.Loss.1,
         Chr_Gene.2 = Other.Arm.Gain.Loss.1,
         Start.2 = Other.Position.Start.1,
         End.2 = Other.Position.End.1,
         Origin.2 = Origin.1,
         Gain_Loss.3 = Other.Gain.Loss.2,
         Chr_Gene.3 = Other.Arm.Gain.Loss.2,
         Start.3 = Other.Position.Start.2,
         End.3 = Other.Position.End.2,
         Origin.3 = Origin.2,
         Gain_Loss.4 = Other.Gain.Loss.3,
         Chr_Gene.4 = Other.Arm.Gain.Loss.3,
         Start.4 = Other.Position.Start.3,
         End.4 = Other.Position.End.3,
         Origin.4 = Origin.3) %>%
# Select only relevant columns
  select(Patient.ID,
         Test.Date,
         Genetic.Status,
         Comments,
         Std.Nomenclature,
         Genomic.Position.Start.Max,
         Genomic.Position.End.Max,
         Test.Method,
         Genome.Browser.Build,
         Result.type,
         Gain_Loss.1,
         Chr_Gene.1,
         Start.1,
         End.1,
         Origin.1,
         Gain_Loss.2,
         Chr_Gene.2,
         Start.2,
         End.2,
         Origin.2,
         Gain_Loss.3,
         Chr_Gene.3,
         Start.3,
         End.3,
         Origin.3,
         Gain_Loss.4,
         Chr_Gene.4,
         Start.4,
         End.4,
         Origin.4,
         Gene,
         Lab,
         Category,
         Test.Verification,
         `Results.Verified.-.Consultant`,
         Karyotype.Start,
         Karyotype.End,
         Array.Version,
         Array.Confirmation.Studies,
         Start.Exon.Intron:Parental.Origin)

# Delete non-informative test results
Genetics %<>%
  filter(Genetic.Status == "Results Verified",
         Test.Method != "Karyotype",
         Test.Method != "Karyotype/FISH",
         Test.Method != "FISH") %>%
# Keep only the latest informative test results
  group_by(Patient.ID) %>%
  filter(Test.Date == last(Test.Date)) %>%
  ungroup

# Small corrections
Genetics$Chr_Gene.2 <- gsub("[pq]$", "", Genetics$Chr_Gene.2)
Genetics$Chr_Gene.3 <- gsub("[pq]$", "", Genetics$Chr_Gene.3)
Genetics$Chr_Gene.4 <- gsub("[pq]$", "", Genetics$Chr_Gene.4)
Genetics$Result.type[Genetics$Test.Method == "Bi-Directional Sequence Analysis" | Genetics$Test.Method == "Sequence Analysis" | Genetics$Test.Method == "Whole exome sequencing"] <- "mutation"

# Select out irrelevant variables
Genetics <- select(Genetics, -Genetic.Status)

Genetics %<>%
  mutate(Test.Date = format(Genetics$Test.Date, format = "%m/%d/%y %I:%M %p"))

write.csv(Genetics, "dataGenetics.csv", na = "", row.names = F)
