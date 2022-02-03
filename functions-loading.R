library(tidyverse)

# Transform an HTML file to a 2 header csv with pipes
html2csv <- function(filename, encoding = "UTF-8")
{
  # Open file as raw UTF-8 text
  file <- file(paste0(filename,".xls"), encoding = encoding)
  html <- readLines(file)
  close(file)

  # Concatenate all lines and preserve newlines
  html <- paste(html,sep = "", collapse = "\n")

  # Keep only the contents of the table
  html <-  sub(".*?<table border=0>(.*?)</table>.*", "\\1",  html)

  # Process th's with colspan attribute (first header row)
  # m is the regexp object matching th's with colspan attribute
  # nb is the number of colspan extracted from the regexp
  while (length(grep("<th.*?colspan=(\\d+).*?>.*?</th>", html)) > 0)
  {
    m <- regexpr("[[:space:]]*<th.*?colspan=(\\d+).*?>[[:space:]]*(.*?)[[:space:]]*</th>[[:space:]]*(?s)", html, perl = T)
    nb <- as.numeric(substr(html, attr(m,"capture.start")[1], attr(m,"capture.start")[1] + attr(m,"capture.length")[1] - 1))
    str <- substr(html, attr(m,"capture.start")[2], attr(m,"capture.start")[2] + attr(m,"capture.length")[2] - 1)
    str <- paste0("|", str, "|", ",")
    html <- sub("[[:space:]]*<th.*?colspan=(\\d+).*?>[[:space:]]*(.*?)[[:space:]]*</th>[[:space:]]*(?s)", paste0(rep(str, nb), collapse = ""), html, perl = T)
  }

  # Process second header row
  html <- gsub("[[:space:]]*</th>[[:space:]]*<th.*?>[[:space:]]*(?s)", '|,|', html, perl = T)
  html <- gsub("[[:space:]]*</?th.*?>[[:space:]]*(?s)",                '|',   html, perl = T)

  # Process body of the table
  html <- gsub("[[:space:]]*</td>[[:space:]]*<td.*?>[[:space:]]*",     '|,|', html)
  html <- gsub("[[:space:]]*</?td.*?>[[:space:]]*",                    '|',   html)

  # Create the correct newlines
  html <- gsub("[[:space:]]*<tr.*?>[[:space:]]*",                      "",    html)
  html <- gsub("[[:space:]]*</tr>[[:space:]]*",                        "\n",  html)

  # Get rid of other tags
  html <- gsub("<.*?>",                                                "",    html)

  # Unescape quotes
  html <- gsub('\\\\+',                                                '',    html)

  # Get rid of trailing commas
  html <- gsub(',\n',                                                  '\n',  html, perl = T)

  # Replace &reg; with ®
  html <- gsub('&reg;',                                                '®',   html)

  # Replace no-break spaces with normal spaces
  html <- gsub(' ',                                                    ' ',   html)

  # Get rid of multiple spaces
  html <- gsub(' +',                                                   ' ',   html)

  # Get rid of trailing whitespaces
  html <- gsub(' +\\|',                                                '|',   html)
  html <- gsub('\\| +',                                                '|',   html)

  # Delete empty fields
  html <- gsub("\\|\\|",                                               "",    html)

  # Write final csv file
  cat(html,file = paste0(filename, ".csv"))
}

# Load a csv file with 2 header rows
read.csv.2header <- function(file, ...)
{
  headers <- read2headers(file)

  # Concatenate and clean the two header lines
  header <- catClean(headers[[1]], headers[[2]])

  # Read the data itself and attribute column names
  data <- read.csv(file, quote = "|", header = F, skip = 2, stringsAsFactors = F, colClasses = "character", ...)
  colnames(data) <- header

  data
}

# Create a template premapping file from a data file
writePremap <- function(datafile, premapfile)
{
  headers <- read2headers(datafile)

  # Concatenate and clean the two header lines
  header <- catClean(headers[[1]], headers[[2]])

  # Re-split the header for the mapping file
  Head1 <- sub("(^.*?)_.*$","\\1", header, perl = T)
  Head2 <- sub("^.*?_(.*$)","\\1", header, perl = T)

  # If only one header, put it in first header
  Head2[Head2 == Head1] <- ""

  # Create all columns for mapping file
  ColNum <- 1:length(Head1)
  premap <- data.frame(ColNum, Head1, Head2, stringsAsFactors = F)
  premap <- mutate(premap, SubFile = "", Evo = "", Reformat = "", VarName = "", Linked = "")
  premap[grepl("\\d+_", premap$Head2), ] <- premap %>%
    filter(grepl("\\d+_", Head2)) %>%
    mutate(Linked = sub("(^\\d+)_.*", "\\1", Head2)) %>%
    mutate(VarName = sub("\\d+_(.*$)", "\\1", Head2))
  premap <- mutate(premap, Header = header)

  write.table(premap, file = premapfile, row.names = F, sep = ",", quote = T, na = "")
}

# Clean and concatenate two headers
catClean <- function(header1, header2)
{
  # Clean variable names
  header1 <- cleanHeader1(header1)
  header2 <- cleanHeader2(header2)

  # Merge the two headers
  header = paste(header1, header2, sep = "_")

  # Clean the merging (trailing "_")
  header <- sub("^_", "", header, perl = T)
  header <- sub("_$", "", header, perl = T)

  # Replace spaces with dots
  header <- gsub(" +", ".", header, perl = T)

  header
}

# Clean the first header
cleanHeader1 <- function(header1)
{
  header1 <- gsub(" +\\(.*?\\)",                                                                                 "", header1, perl = T)
  header1 <- gsub("Please (enter|select) either pounds( and ounces)? or kilograms\\.",                           "", header1, perl = T)
  header1 <- gsub("Please enter either feet and inches or centimeters\\.",                                       "", header1, perl = T)
  header1 <- gsub("Please answer the following questions\\.",                                                    "", header1, perl = T)
  header1 <- gsub("If you answer yes to any of the following questions, please select the age of occurrence\\.", "", header1, perl = T)

  header1
}

# Clean the second header
cleanHeader2 <- function(header2)
{
  header2 <- gsub("^Responses$",                          "",        header2, perl = T)
  header2 <- gsub(" - APGAR score$",                      "",        header2, perl = T)
  header2 <- gsub(" - Frequency$",                        "",        header2, perl = T)
  header2 <- gsub("([^\\d])\\.Age( at milestones?)?$",    "\\1",     header2, perl = T)
  header2 <- gsub("\\?[[:space:]]*\\.",                   "?_",      header2, perl = T)
  header2 <- gsub("\"",                                   "",        header2, perl = T)
  header2 <- gsub("^OpenText$",                           "",        header2, perl = T)
  header2 <- gsub("(\\(during the day\\))\\.",            "\\1_",    header2, perl = T)
  header2 <- gsub("(\\(out of any type of vessel\\))\\.", "\\1_",    header2, perl = T)
  header2 <- gsub("(\\d)\\.(\\w)",                        "\\1_\\2", header2, perl = T)
  header2 <- gsub("\\.(Symptoms|Age at Symptoms|Status)","_\\1",     header2, perl = T)

  header2
}

# Read headers from a two headers file and propagate the first one to all corresponding columns
read2headers <- function(file)
{
  header1 <- scan(file, what = character(), nlines = 1, sep = ",", quote = "|")
  header2 <- scan(file, what = character(), nlines = 1, sep = ",", quote = "|", skip = 1)

  # Return the two headers as a list
  headers <- list(header1, header2)

  headers
}

