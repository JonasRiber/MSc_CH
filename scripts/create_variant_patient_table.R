




### author: Jonas Riber Jørgensen
### date: 24-04-2025
### desc: combine the CH+ annotated variants into a table with var in rows and patient in col


library(tidyverse)
library(data.table)


setwd("/Path/to/working/directory")
if(!interactive()) pdf(NULL) # suppress Rplots.pdf file

### setup args
args <- commandArgs(trailingOnly = TRUE)

#check if a file path was provided
if (length(args) < 2) {
  stop("Usage: Rscript mut_signature.R <input_directory.tsv> <output_directory>") # input: takes all files in the directory that ends with _full.tsv
}

folder_path <- args[1]   #First argument is the file path
output_dir <- args[2]  #Second argument is output directory path


#check if in and out location exists
if (!file.exists(folder_path)) {
  stop(paste("File not found:", folder_path))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}







#list of all .tsv files
tsv_files <- list.files(folder_path, pattern = "\\marked_filtered.tsv$", full.names = TRUE)

nr_files <- length(tsv_files)



#function to read each file and add a column with the id from filename
read_and_tag <- function(file) {
  read.delim(file, stringsAsFactors = FALSE) %>%
    mutate(id = strsplit(basename(file), "_")[[1]][1])
}


empty_patients <- character(0)

read_and_tag <- function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)
  patient_id <- strsplit(basename(file), "_")[[1]][1]
  
  if (nrow(df) == 0) {
    message("Empty file detected, creating dummy structure for: ", basename(file))
    
    # add to list
    empty_patients <<- c(empty_patients, patient_id) 
    
    # Manually create an empty data frame with correct column types
    df <- data.frame(
      Chr = character(0),
      Start = integer(0),
      End = integer(0),
      Ref = character(0),
      Alt = character(0),
      AA_change = character(0),
      Func = character(0),
      Gene = character(0),
      Alt_count = numeric(0),
      stringsAsFactors = FALSE
    )
    
    
  }
  
  # Add the patient ID (safe even for 0-row frames)
  df <- df %>% mutate(id = patient_id)
  
  # Ensure Alt column is character if present (just in case)
  if ("Alt" %in% names(df)) {
    df$Alt <- as.character(df$Alt)
  }
  
  return(df)
}





#apply the function to all files and combine them
data <- bind_rows(lapply(tsv_files, read_and_tag))


data <- data |> filter(Alt_count<Somatic_threshold)



# Get all unqiue variants
unique_variants <- data |> 
  distinct(Chr, Start, End, Ref, Alt, AA_change, Func, Gene)

unique_patients <- data |> distinct(id)

all_combinations <- tidyr::crossing(unique_variants, unique_patients)
all_combinations$Start <- all_combinations$Start-1 # need to use depth caller again, uses 0-indexing

patients_data <- split(all_combinations, all_combinations$id)

dir.create(output_dir, showWarnings = FALSE) 

#write 0-position files
lapply(names(patients_data), function(patient) {
  write.table(patients_data[[patient]], paste0(output_dir, patient, "_0_positions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
})



### handle empty files
#use last combination file as template
template_file <- file.path(output_dir, tail(list.files(output_dir, pattern = "_0_positions.tsv$"), 1))



# create missing output files as copies of the last good one
for (pid in empty_patients) {
  new_file <- file.path(output_dir, paste0(pid, "_0_positions.tsv"))
  file.copy(template_file, new_file)
  message("Created placeholder file for patient: ", pid)
}




