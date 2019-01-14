require("dplyr")

processRanges <- function(genetics)
{
  hg38 <- read.delim("refGene.txt.hg38", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  # Prepare the data frame to contain only coordinates by transforming gene information and mutation information
  genetics <- select(genetics, Patient.ID, Genome.Browser.Build, Result.type, Gain_Loss, Chr_Gene, Start, End)
  genetics$Start <- as.numeric(genetics$Start)
  genetics$End <- as.numeric(genetics$End)
  for (i in 1:nrow(genetics))
  {
    if (genetics$Result.type[i] == "gene")
    {
      genetics$Genome.Browser.Build[i] <- "GRCh38/hg38"
      genetics$Start[i]                <-     as.numeric(hg38$txStart[hg38$name2 == genetics$Chr_Gene[i]][1])
      genetics$End[i]                  <-     as.numeric(hg38$txEnd  [hg38$name2 == genetics$Chr_Gene[i]][1])
      genetics$Chr_Gene[i]             <- sub("chr", "", hg38$chrom  [hg38$name2 == genetics$Chr_Gene[i]][1])
    }
    else if (genetics$Result.type[i] == "mutation")
    {
      genetics$Genome.Browser.Build[i] <- "GRCh38/hg38"
      genetics$Start[i]                <- genetics$Start[i] + as.numeric(hg38$txStart[hg38$name2 == genetics$Chr_Gene[i]][1])
      genetics$End[i]                  <- genetics$End[i]   + as.numeric(hg38$txStart[hg38$name2 == genetics$Chr_Gene[i]][1])
      genetics$Chr_Gene[i]             <- sub("chr", "",      hg38$chrom  [hg38$name2 == genetics$Chr_Gene[i]][1])
    }
  }
  
  liftOver(genetics)
}

liftOver <- function(genetics)
{
  # Create input bed files
  genetics <- add_rownames(genetics)
  genetics$rowname <- as.numeric(genetics$rowname)
  genetics %>%
    filter(Genome.Browser.Build == "NCBI35/hg17") %>%
    select(Chr_Gene, Start, End, rowname) %>%
    mutate(Chr_Gene = paste0("chr", Chr_Gene)) -> bed17

  genetics %>%
    filter(Genome.Browser.Build == "NCBI36/hg18") %>%
    select(Chr_Gene, Start, End, rowname) %>%
    mutate(Chr_Gene = paste0("chr", Chr_Gene)) -> bed18

  genetics %>%
    filter(Genome.Browser.Build == "NCBI37/hg19" | Genome.Browser.Build == "GRCh37/hg19") %>%
    select(Chr_Gene, Start, End, rowname) %>%
    mutate(Chr_Gene = paste0("chr", Chr_Gene)) -> bed19

  write.table(bed17, file = "bed17", sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(bed18, file = "bed18", sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(bed19, file = "bed19", sep = "\t", row.names = F, col.names = F, quote = F)

  # Run liftOver
  system("./liftOver bed17 hg17ToHg19.over.chain out17 unmap17 -multiple")
  system("./liftOver out17 hg19ToHg38.over.chain out17b unmap17b -multiple")
  system("./liftOver bed18 hg18ToHg38.over.chain out18 unmap18 -multiple")
  system("./liftOver bed19 hg19ToHg38.over.chain out19 unmap19 -multiple")

  out17 <- read.delim("out17b", header = F, stringsAsFactors = F)
  out18 <- read.delim("out18",  header = F, stringsAsFactors = F)
  out19 <- read.delim("out19",  header = F, stringsAsFactors = F)

  # Read the liftOver results and keep only mapped regions on the original chromosome
  out17 <- filter(out17, V1 == paste0("chr",genetics$Chr_Gene[V4]))
  out18 <- filter(out18, V1 == paste0("chr",genetics$Chr_Gene[V4]))
  out19 <- filter(out19, V1 == paste0("chr",genetics$Chr_Gene[V4]))

  # Join together the overlapping regions on each chromosome
  data <- matrix(ncol = 4, nrow = 0)
  for (pat in unique(out17$V4))
  {
    depth <- 0
    start <- 0
    end   <- 0

    pos <- sort(unlist(out17[out17$V4 == pat, 2:3]))
    for (posi in 1:length(pos))
    {
      if (grepl("^V2", names(pos[posi])))
      {
        depth <- depth + 1
        if (depth == 1)
          start <- pos[posi]
      }
      else
        depth <- depth - 1

      if (depth == 0)
      {
        end <- pos[posi]
        data <- rbind(data, c(pat, sub("chr", "", unique(out17$V1[out17$V4 == pat])), start, end))
      }
    }
  }
  for (pat in unique(out18$V4))
  {
    depth <- 0
    start <- 0
    end   <- 0

    pos <- sort(unlist(out18[out18$V4 == pat, 2:3]))
    for (posi in 1:length(pos))
    {
      if (grepl("^V2", names(pos[posi])))
      {
        depth <- depth + 1
        if (depth == 1)
          start <- pos[posi]
      }
      else
        depth <- depth - 1

      if (depth == 0)
      {
        end <- pos[posi]
        data <- rbind(data, c(pat, sub("chr", "", unique(out18$V1[out18$V4 == pat])), start, end))
      }
    }
  }
  for (pat in unique(out19$V4))
  {
    depth <- 0
    start <- 0
    end   <- 0

    pos <- sort(unlist(out19[out19$V4 == pat, 2:3]))
    for (posi in 1:length(pos))
    {
      if (grepl("^V2", names(pos[posi])))
      {
        depth <- depth + 1
        if (depth == 1)
          start <- pos[posi]
      }
      else
        depth <- depth - 1

      if (depth == 0)
      {
        end <- pos[posi]
        data <- rbind(data, c(pat, sub("chr", "", unique(out19$V1[out19$V4 == pat])), start, end))
      }
    }
  }

  data <- data.frame(data, stringsAsFactors = F)
  names(data) <- c("rowname", "Chr_Gene", "Start", "End")
  data$rowname <- as.numeric(data$rowname)
  data$Genome.Browser.Build <- "GRCh38/hg38"
  data$Result.type <- "coordinates"
  data$Gain_Loss <- genetics$Gain_Loss[data$rowname]
  data$Patient.ID <- genetics$Patient.ID[data$rowname]

  data <- select(data, rowname, Patient.ID, Genome.Browser.Build, Result.type, Gain_Loss, Chr_Gene, Start, End)

  genetics <- filter(genetics, !(rowname %in% data$rowname))
  genetics <- rbind(genetics, data)
  genetics <- arrange(genetics, rowname)

  # unlink(c("bed*","out*","unmap*"))

  select(genetics, -rowname)
}

downloadExternalFiles <- function()
{
  # LiftOver
  system("wget -O liftOver http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O - | gunzip > hg19ToHg38.over.chain")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz -O - | gunzip > hg18ToHg38.over.chain")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz -O - | gunzip > hg17ToHg19.over.chain")
}
