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

html2csv("dataClinical")
html2csv("dataDevelopmental")
html2csv("dataAdult")
# html2csv("dataGenetic")

library(tidyverse)
library(readxl)

read_xlsx("dataGenetic.xlsx") %>%
  write_csv("dataGenetic.csv")
