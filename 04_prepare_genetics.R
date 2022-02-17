library(tidyverse)

# Read genetic file and select relevant columns, in correct ordering
# Discard columns that could contain identifiable data
# Discard columns that don't present analytical usability
read_csv("dataGenetic.csv") -> genetic

genetic %>%
  select(Patient.ID,
         GeneticStatus,
         SHANK3Involved,
         PathogenicSHANK3Defect,
         # CategoryChr22,
         # SubCategoryChr22,
         Mosaicism,
         Origin,
         TestVerification,
         # ParentTestMethod,
         # ParentResults,
         # starts_with("Karyotype"),
         # starts_with("FISH"),
         Array,
         starts_with("Array1"),
         # starts_with("Sequencing"),
         Sequencing,
         SequencingGenomicChange =  GenomicChange,
         SequencingCodingDNAChange = CodingDNAChange,
         SequencingProteinChange = ProteinChange,
         SequencingEffect = Effect,
         SequencingClinicalSignificance = ClinicalSignificance,
         # SequencingOtherAbnormalities,
         -ends_with("Lab"),
         -ends_with("Version"),
         -ends_with("Method"),
         -contains("Other"),
         -contains("MAX")) %>%

  # Fix some category levels.
  # (-> Mosaicism and Origin should be fixed upstream)
  mutate(Mosaicism        = Mosaicism        %>% fct_collapse("Yes"                   = c("Yes", "Yes (duplication)")),
         Origin           = Origin           %>% fct_collapse("De novo"               = c("De novo", "De Novo"),
                                                              "Unknown"                = c("Unknown", "Unknown (see CommentsCB)"),
                                                              "Unknown - Not maternal" = c("Unknown - Not maternal", "Unknown (Not maternal)"),
                                                              "Unknown - Not paternal" = c("Unknown - Not paternal", "Unknown - Not Paternal", "Unknown (Not paternal)", "Unknown (not paternal)")),
         TestVerification = TestVerification %>% fct_collapse("Complete"              = genetic$TestVerification %>% keep(~str_detect(., "^Complete")) %>% unique,
                                                              "Incomplete"             = genetic$TestVerification %>% keep(~str_detect(., "^Incomplete")) %>% unique)) %>%
# Ensure years are numeric
  mutate_at(vars(ends_with("Year")), as.numeric) %>%
# Keep only verified results (should select all results now)
  filter(GeneticStatus == "Results verified" | GeneticStatus == "Results Verified")  -> genetic

# #Prepare names for reshaping and transmart mapping
# genetic %>%
#   set_names(function(x)
#             {
#               x %>%
#                 str_replace("^Array1_", "Array.") %>%
#                 str_replace("Chr22Abnl_", "") %>%
#                 str_replace("MIN", "") %>%
#                 str_replace("(\\d)_GRCh38hg38", ".latest\\1") %>%
#                 str_replace("(\\d)$", "_\\1")
#             }) -> genetic

# Prepare names for reshaping and transmart mapping
rename <- function(x)
  {
    x %>%
      str_replace("^Array1_", "Array.") %>%
      str_replace("Chr22Abnl_", "") %>%
      str_replace("MIN", "") %>%
      str_replace("(\\d)_GRCh38hg38", ".latest\\1") %>%
      str_replace("(\\d)$", "_\\1")
  }
genetic_name <- rename(colnames(genetic))
genetic %>% set_names(genetic_name) -> genetic

# Reshape
genetic %>%
  filter(Array == "Yes") %>%
  reshape(direction = "long",
          varying = genetic %>% names %>% str_detect("_\\d$"),
          sep = "_",
          timevar = "Array.num") %>%
  filter(!Array.Type %>% is.na) %>%
  mutate(Result.type = "Array") %>%
  select(-id, -starts_with("Sequencing"), -Array) %>%
  bind_rows(genetic %>%
              filter(Sequencing == "Yes") %>%
              mutate(Result.type = "Sequencing") %>%
              select(-starts_with("Array"), -Sequencing)) -> genetic

genetic %>%
  write_csv("dataGenetic.csv")

