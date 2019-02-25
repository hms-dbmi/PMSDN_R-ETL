README
======

R scripts to load and clean PMSIR raw data files.

This repo uses `packrat` for R libraries dependencies management to ensure reproducibility of execution. After cloning, launching R within the repo should auto-start the installation of R libraries dependencies in the **packrat** subdir of the repo.  
**If this step does not start automatically,** run `packrat::restore()` to manually download all dependencies.

Step 0
------

Get the files from PMSIR.  
Four input files are obtained from the registry:

* Clinical Questionnaire
* Developmental Questionnaire
* Adolescent and Adults Questionnaire
* Genetic test results

The three first files are HTML tables in .xls files. **Do not open these files with excel !**

They must first be renamed as following:
* dataClinical.xls for the Clinical Questionnaire
* dataDevelopmental.xls for the Developmental Questionnaire
* dataAdult.xls for the Adolescent and Adults Questionnaire
* dataGenetic.xlxs for the Genetic test results

and placed in the scripts folder.

Optionally run *00_opt_download_files.sh* to download the necessary files to run `liftOver`:
  + the `liftOver` tool itself
  + hg17ToHg19.over.chain.gz : mapping file from Human Genome Assembly Hg17 to Hg19
  + hg18ToHg38.over.chain.gz : mapping file from Human Genome Assembly Hg18 to Hg38
  + hg19ToHg38.over.chain.gz : mapping file from Human Genome Assembly Hg19 to Hg38

The `liftOver` tool is used to convert chromosome coordinates to the latest Human Genome Assembly (GRCh38/hg38):  
* hg17 -> hg19 -> hg38
* hg18 -> hg38
Some coordinates still cannot be converted, so the result is a mix of different human genome assemblies, with the fewest not up to date as possible.

Step 1
------

The files must be converted to the csv format.  
This is done using the first script: *01_prepare_files.sh*  
This script calls the *01_prepare.R* R script.  
Four new files are created:
* dataClinical.csv
* dataDevelopmental.csv
* dataAdult.csv
* dataGenetic.csv

Please read the source file as it contains specific information about the current raw data files we obtained.

All but the dataGenetic file are **UTF-8 encoded**, **comma (,) separated**, with **pipes (|) as quotes** to delimit text fields.  
You can check the csv files by opening them in LibreOffice (excel doesn't accept pipes as text delimiters)  
This format keeps the double header format from the registry, which can be used to regenerate premapping files in a later optional step.

![](Docs/libreoffice.png)

Step 2 (optional)
-----------------

Regenerate the premapping files.  
The *02_opt_create_premap_files.sh* script uses the data files to generate the empty premapping files.  
This script calls the *02_opt_create_premap_files.R* R script, which in turn uses functions from the *functions-loading.R* script.  
The script creates three premapping files, one for each data file:
* premapClinical.csv
* premapDevelopmental.csv
* premapAdult.csv

The premapping files are used to tell how to process each variable and where to put it in the ontology.  
Already filled-in premapping files are on the repo and will work as long as the structure of the registry doesn't change.  
In case there are modifications to make to the processing, these premapping files can be amended.  
In case the registry changes (new variables, renamed variables, etc.), the premapping files can be generated again from scratch.
They then have to be completed again, mainly using the previous ones as a template.  
A more in-depth manual on how to use the premapping files and add processing code to the scripts can be found in a dedicated section of this readme.

Step 3
------

Some entries are duplicated, or test accounts.  
This step deletes true duplicated lines, lines from the test accounts, and merges the identity of the duplicated accounts.
It then saves the files as regular CSV files with a single header row.

The duplicated and tests accounts must be given as a `duplicates_and_tests.csv` csv file, with the following structure:

| type | id1                         | id2                          |
|------|-----------------------------|------------------------------|
| test | _id of the test account_    |                              |
| dup  | _id of the first duplicate_ | _id of the second duplicate_ |

Some caveat apply when choosing which id to put in the **id1** and **id2** columns. **id2** should be the id in use in the genetics file, and/or the one with the most recent phenotypic data.

Step 4
------

Prepare the genetics file.  
This step applies specific preparation for the genetics file:

- selecting relevant variables
- excluding variables containing potentially identifiable data
- pre-cleaning variables (factor levels management, etc.)
- reshaping the data for use in analysis and integration in transmart (one row per unique genetic results, multiple rows for each patient with multiple results)

Genetic data at this point is still in the original presentation.


Step 5
------

Each file is processed, generating outputs for analysis and for integration into transmart.
The *05_clean_data.sh* scripts calls the *05_clean_data.R* R script, which in turn uses functions from:
* *functions-loading.R* to load the data and premapping files
* *functions-process.R* to process each level of depth of the data
* *functions-mapping.R* to generate the Kettle mapping file
* *functions-reformatting.R* to reformat and clean variables according to the rules in the premapping files
* *functions-genes.R* to process genetic data (data files and mapping file)

This final script creates two directories:
- `output_human/` directory containing data files usable for analysis
- `output_transmart/` directory containing the cleaned data files (one for each section of each questionnaire) and the Kettle mapping file for integration into transmart.

The files **[Cleaning data overview.pdf](Docs/Cleaning data overview.pdf)** and **[PMSIR Data Management Algorithm Detailed.pdf](Docs/PMSIR Data Management Algorithm Detailed.pdf)** illustrate how this step works internally.  
The source files (*.svg) for these diagrams are provided.
***

Premapping files
================

The premapping files define:
* what are the useful variables in the data files
* where they are located
* where they should go in the ontology tree
* if the variable is an historic variable or an evolutive one
* if and what reformatting is needed
* how different variables are linked together, which is used by the reformatting option

The structure of the ontology tree in i2b2/tranSMART is as follows:
* Questionnaire
  + SubFile
    + First header
      + Second header

For example:
* Clinical
  + Ears - Hearing
    + Has the patient had any of the following hearing tests?
      + Behavioral audiometry
      + Tympanogram
    + Has the patient had ear tubes
  + Mouth - Dental
    + Has the patient ever ground his or her teeth?
      + Yes, during the day
      + Yes, during the night
      + No
      + Unsure

The structure of the file is as follows:

ColNum | Head1 | Head2 | SubFile | Evo | Reformat | VarName | Linked | Header
------ | ----- | ----- | ------- | --- | -------- | ------- | ------ | ------
The column number in which to find the variable in the data file | The first header | The second header | The "subfile" in which it should reside | The evolutive status | The reformatting function to use | The new variable name after reformatting | A number to show links between variables | The complete name of the variable found in the data files
Auto   | Auto  | Auto  | Manual | Manual | Manual | Auto    | Auto   | Auto
 |       |       | From the registry structure, absent from the export files | From epidemiologist expertise. Any variable refering to the current status of the patient | Function to use to reformat the variable | | |
 |       |       | Plain text string | 1 if evolutive, nothing if not | name of the R function to use | | |

If you want to add new reformatting options, write an R function in *functions-reformatting.R* with a unique name that will be used to reformat the data, such as the *refactor* function.
A prototype for such a function (here it does nothing and just uses the data as-is) is:
```R
new_reformat_fn <- function(data, premap)
{
  # The 'premap' object is the subset of the premapping corresponding to the
  # current First Header level of the ontology processed
  # The 'data' object is the subset of data corresponding to the subset of the premap object passed
  # So in theory, only one reformatting function can be used for one "First Header" subset of variables.
  # Which cells are marked with the reformatting option, and the content of the
  # Linked cell can be used as hints to control the behavior of the function
  # The whole function could be written as just "data" to return the data
  # exactly as-is. But here I show how to use relevant information from the two passed objects

  # Create new data frame to contain transformed/curated data
  # This object must always contain the three variables Patient.ID, Survey.Time and Birthdate
  data2 <- select(data, Patient.ID, Survey.Time, Birthdate)

  # Create the variables in the output object, here with the same names as
  # originally, which can be found in the premap object
  # This is where variable selection and/or creation can be done.
  varnames <- premap$Header 

  # Here we just copy the contents from the original data to the output data2.
  # This is where transformation can be done.
  data2[varnames] <- data[varnames]

  # Return the new object
  data2
}
```
***
