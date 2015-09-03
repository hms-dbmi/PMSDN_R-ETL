require("dplyr")

processRanges <- function(genetics)
{
  hg38 <- read.delim("refGene.txt.hg38", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  # Prepare the data frame to contain only coordinates by transforming gene information and mutation information
  genetics <- select(genetics, Patient.ID, Genome.Browser.Build, Result.type, Gain_Loss, Chr_Gene, Start, End)
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

processGenes <- function(genetics, bin = F)
{
  # Get a list of all involved genes
  genes <- getGenes(genetics)
  
  # Create the data frame to hold the annotated genetic data
  Genetics_genes <- data_frame(Patient.ID = unique(genetics$Patient.ID))
  Demographics <- processDemographics(noOutput = T)
  Demographics$Patient.ID <- as.numeric(Demographics$Patient.ID)
  Genetics_genes <- left_join(Genetics_genes, Demographics[c("Patient.ID", "Gender")])
  for (gene in genes$name)
  {
    if (bin)
    {
      if (genes$chrom[genes$name == gene] == "Y")
      {
        Genetics_genes[[gene]][Genetics_genes$Gender == "Male"]   <- F
        Genetics_genes[[gene]][Genetics_genes$Gender == "Female"] <- NA
      }
      else
        Genetics_genes[[gene]] <- F
    }
    else
    {
      if (genes$chrom[genes$name == gene] == "Y")
      {
        Genetics_genes[[gene]][Genetics_genes$Gender == "Male"]   <- 1
        Genetics_genes[[gene]][Genetics_genes$Gender == "Female"] <- 0
      }
      else if (genes$chrom[genes$name == gene] == "X")
      {
        Genetics_genes[[gene]][Genetics_genes$Gender == "Male"]   <- 1
        Genetics_genes[[gene]][Genetics_genes$Gender == "Female"] <- 2
      }
      else
        Genetics_genes[[gene]] <- 2
    }
  }
  
  # Extract the information from the raw genetic test reports into the data frame
  extractGenes(genetics, Genetics_genes, bin)
}

processPathways <- function(genetics, genetics_genes)
{
  # Enrich genes with pathways annotation
  genes <- getGenes(genetics)
  genes <- getPathways(genes)
  
  # Create the data frame to hold the annotated genetic data
  Genetics_pathways <- data.frame(Patient.ID = unique(genetics$Patient.ID))
  Genetics_pathways[unique(sort(genes$pathway))] <- 0
  
  extractPathways(genetics_genes, genes, Genetics_pathways)
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

  unlink(c("bed*","out*","unmap*"))

  select(genetics, -rowname)
}

getGenes <- function(genetics)
{
  # Read refGene position tables for all genome assemblies
  hg17 <- read.delim("refGene.txt.hg17", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg18 <- read.delim("refGene.txt.hg18", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg19 <- read.delim("refGene.txt.hg19", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg38 <- read.delim("refGene.txt.hg38", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))

  genes <- data.frame(name = character(0), chrom = character(0))

  for (row in 1:nrow(genetics))
  {
    if (genetics$Result.type[row] == "mutation" | genetics$Result.type[row] == "gene")
      genes <- rbind(genes, data.frame(name = genetics$Chr_Gene[row], chrom = "", stringsAsFactors = F))
    else if (genetics$Result.type[row] == "coordinates")
    {
      if      (genetics$Genome.Browser.Build[row] == "GRCh37/hg19")
        genome <- hg19
      else if (genetics$Genome.Browser.Build[row] == "GRCh38/hg38")
        genome <- hg38
      else if (genetics$Genome.Browser.Build[row] == "NCBI35/hg17")
        genome <- hg17
      else if (genetics$Genome.Browser.Build[row] == "NCBI36/hg18")
        genome <- hg18
      else
        genome <- hg18

      if (genetics$Gain_Loss[row] == "Loss")
      {
        name <- genome$name2[((as.numeric(genome$txEnd) > genetics$Start[row] & as.numeric(genome$txEnd) < genetics$End[row]) | (as.numeric(genome$txStart) > genetics$Start[row] & as.numeric(genome$txStart) < genetics$End[row]) | (as.numeric(genome$txStart) < genetics$Start[row] & as.numeric(genome$txEnd) > genetics$End[row])) & genome$chrom == paste0("chr", genetics$Chr_Gene[row])]
        if (length(name) > 0)
          genes <- rbind(genes, data.frame(name, chrom = genetics$Chr_Gene[row], stringsAsFactors = F))
      }
      else
      {
        name <- genome$name2[as.numeric(genome$txStart) > genetics$Start[row] & as.numeric(genome$txEnd) < genetics$End[row] & genome$chrom == paste0("chr", genetics$Chr_Gene[row])]
        if (length(name) > 0)
          genes <- rbind(genes, data.frame(name, chrom = genetics$Chr_Gene[row], stringsAsFactors = F))
      }

    }
  }

  genes <- genes %>% arrange(name, chrom) %>% distinct
  genes$chrom[genes$chrom == ""] <- sub("chr", "", (hg38[hg38$name2 %in% genes$name[genes$chrom == ""], c("chrom","name2")] %>% arrange(name2, chrom) %>% distinct)$chrom)
  genes <- genes %>% arrange(name, chrom) %>% distinct

  genes
}

extractGenes <- function(genetics_pre, genetics_post, bin)
{
  # Read refGene position tables for all genome assemblies
  hg17 <- read.delim("refGene.txt.hg17", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg18 <- read.delim("refGene.txt.hg18", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg19 <- read.delim("refGene.txt.hg19", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))
  hg38 <- read.delim("refGene.txt.hg38", stringsAsFactors = F, header = F, col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"))

  for (row in 1:nrow(genetics_pre))
  {
    patient <- genetics_pre$Patient.ID[row]
    Chr_Gene <- genetics_pre$Chr_Gene[row]

    if (genetics_pre$Result.type[row] == "mutation")
    {
      if (bin)
        genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] <- T
      else
        genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] <- genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] - 1
    }
    else if (genetics_pre$Result.type[row] == "gene")
    {
      if (bin)
        genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] <- T
      else
      {
        if (genetics_pre$Gain_Loss[row] == "Gain")
          genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] <- genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] + 1
        else
          genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] <- genetics_post[genetics_post$Patient.ID == patient, Chr_Gene] - 1
      }
    }
    else if (genetics_pre$Result.type[row] == "coordinates")
    {
      if      (genetics_pre$Genome.Browser.Build[row] == "GRCh37/hg19")
        genome <- hg19
      else if (genetics_pre$Genome.Browser.Build[row] == "GRCh38/hg38")
        genome <- hg38
      else if (genetics_pre$Genome.Browser.Build[row] == "NCBI35/hg17")
        genome <- hg17
      else if (genetics_pre$Genome.Browser.Build[row] == "NCBI36/hg18")
        genome <- hg18
      else
        genome <- hg18

      if (genetics_pre$Gain_Loss[row] == "Loss")
      {
        genes <- unique(genome$name2[((as.numeric(genome$txEnd) > genetics_pre$Start[row] & as.numeric(genome$txEnd) < genetics_pre$End[row]) | (as.numeric(genome$txStart) > genetics_pre$Start[row] & as.numeric(genome$txStart) < genetics_pre$End[row]) | (as.numeric(genome$txStart) < genetics_pre$Start[row] & as.numeric(genome$txEnd) > genetics_pre$End[row])) & genome$chrom == paste0("chr", Chr_Gene)])
        if (bin)
          genetics_post[genetics_post$Patient.ID == patient, genes] <- T
        else
          genetics_post[genetics_post$Patient.ID == patient, genes] <- genetics_post[genetics_post$Patient.ID == patient, genes] - 1
      }
      else
      {
        genes <- unique(genome$name2[as.numeric(genome$txStart) > genetics_pre$Start[row] & as.numeric(genome$txEnd) < genetics_pre$End[row] & genome$chrom == paste0("chr", Chr_Gene)])
        if (bin)
          genetics_post[genetics_post$Patient.ID == patient, genes] <- T
        else
          genetics_post[genetics_post$Patient.ID == patient, genes] <- genetics_post[genetics_post$Patient.ID == patient, genes] + 1
      }
    }
  }

  genetics_post
}

getPathways <- function(genes)
{
  kegg_genes <- read.delim("KEGG_genes.txt", header = F, stringsAsFactors = F)
  kegg_pathways <- read.delim("KEGG_pathways.txt", header = F, stringsAsFactors = F)
  kegg_links <- read.delim("KEGG_link_genes_pathways.txt", header = F, stringsAsFactors = F)

  genes$kegg_gene <- NA
  for (gene in genes$name)
  {
    kegg <- kegg_genes$V1[grep(paste0("^", gene, "[,;]"), kegg_genes$V2)]
    genes$kegg_gene[genes$name == gene] <- ifelse(length(kegg) > 0, kegg, NA)
  }

  genes <- left_join(genes, kegg_links[kegg_links$V1 %in% genes$kegg_gene,], by = c("kegg_gene" = "V1"))
  genes <- rename(genes, kegg_pathway = V2)
  genes <- left_join(genes, kegg_pathways, by = c("kegg_pathway" = "V1"))
  genes <- rename(genes, pathway = V2)
  genes$pathway <- sub(" - Homo sapiens \\(human\\)$","",genes$pathway)

  genes %>% filter(!is.na(pathway))
}

extractPathways <- function(genetics_genes, genes, genetics_post)
{
  for (patient in genetics_genes$Patient.ID)
  {
    for (gene in unique(genes$name))
    {
      # Getting expected number of copies depending on gender and chromosome
      if ((genetics_genes$Gender == "Male" && genes$chrom[genes$name == gene] == "X") | (genetics_genes$Gender == "Male" && genes$chrom[genes$name == gene] == "Y"))
        normal <- 1
      else
        normal <- 2

      if (is.na(genetics_genes[[gene]][genetics_genes$Patient.ID == patient]))
        NA
      else if (genetics_genes[[gene]][genetics_genes$Patient.ID == patient] != normal)
        genetics_post[genetics_post$Patient.ID == patient, genes$pathway[genes$name == gene]] <- genetics_post[genetics_post$Patient.ID == patient, genes$pathway[genes$name == gene]] + 1
    }
  }

  genetics_post
}

downloadExternalFiles <- function()
{
  # KEGG
  system("wget -O KEGG_genes.txt http://rest.kegg.jp/list/hsa")
  system("wget -O KEGG_pathways.txt http://rest.kegg.jp/list/pathway/hsa")
  system("wget -O KEGG_link_genes_pathways.txt http://rest.kegg.jp/link/pathway/hsa")

  # LiftOver
  system("wget -O liftOver http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O - | gunzip > hg19ToHg38.over.chain")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz -O - | gunzip > hg18ToHg38.over.chain")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz -O - | gunzip > hg17ToHg19.over.chain")

  # RefGene
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz -O - | gunzip > refGene.txt.hg38")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz -O - | gunzip > refGene.txt.hg19")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz -O - | gunzip > refGene.txt.hg18")
  system("wget http://hgdownload.soe.ucsc.edu/goldenPath/hg17/database/refGene.txt.gz -O - | gunzip > refGene.txt.hg17")
}
