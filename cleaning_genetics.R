source("functions-process.R")

# Read the curated genetic data
Genetics <- read.csv("dataGenetics.csv", stringsAsFactors = F, na.strings = "")

# Reshape the data frame to prepare it for the coordinates conversion
varying = grep("\\d$", names(Genetics), value = T)
Genetics <- reshape(Genetics, direction = "long", varying = varying, sep = ".", idvar = "Patient.ID", timevar = "Test.nb") %>%
            arrange(Patient.ID, Test.nb) %>%
            filter(!is.na(Gain_Loss))

# ==== Ranges ====
##push("Chromosomal alterations")
Genetics_ranges   <- processRanges(Genetics)

# Keep only GRCh38/hg38 deletion ranges
Genetics_ranges <- filter(Genetics_ranges, Genome.Browser.Build == "GRCh38/hg38") %>% select(-Genome.Browser.Build)
Genetics_ranges <- group_by(Genetics_ranges, Patient.ID, Chr_Gene) %>% mutate(N = row_number()) %>% ungroup() %>% ungroup()

data2 <- distinct(Genetics_ranges[c("Patient.ID", "Result.type")])
for (chr in unique(Genetics_ranges$Chr_Gene))
{
  ranges <- filter(Genetics_ranges, Chr_Gene == chr) %>% select(-Chr_Gene, -Result.type)
  for (n in unique(ranges$N))
  {
    ranges2 <- filter(ranges, N == n)
    names(ranges2) <- c("Patient.ID", paste0("Chr", chr, "_", n,"_",c("Gain_Loss","Start","End")), "N")
    data2 <- merge(data2, select(ranges2, -N), by = "Patient.ID", all = T)
  }
}

##pop()



# ==== Genes ====
Genetics_genes    <- processGenes(Genetics)

# ==== Pathways ====
Genetics_pathways <- processPathways(Genetics, Genetics_genes)
