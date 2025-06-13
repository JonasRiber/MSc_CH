

#### Author: Jonas Riber Jørgensen
#### date: 14/02/2025
#### description: R-script to process and visualize snpcalls

### initialization
library(tidyverse)
library(data.table)

setwd("/path/to/working/directory")

if(!interactive()) pdf(NULL) # suppress Rplots.pdf file


### setup args
args <- commandArgs(trailingOnly = TRUE)

#check if a file path was provided
if (length(args) < 2) {
  stop("Usage: Rscript filtercalls.R <input_file.tsv> <output_path>")
}

file_path <- args[1]   #First argument is the file path
output_path <- args[2]  #Second argument is output directory path


#check if input location exists
if (!file.exists(file_path)) {
  stop(paste("File not found:", file_path))
}



#read input
data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#colnames(data)

setnames(data, old = c("X.CHROM", "START","END", "REF", "ALT", "MTYPE", "N_CON", "N_ALT", "N_ALL", 
                       "VAF", "X"), 
         new = c("chr", "start", "end", "ref", "alt", "Mtype", "n_con", "n_alt", "depth", 
                 "VAF", "somatic_threshold"))
setDT(data)
#head(data)


file_name <- tools::file_path_sans_ext(basename(file_path))
id <- print(str_extract(string = file_name, pattern = "ID_PATTERN"))



log_file <- paste0("data/",id, "/", id, "_filter_log.txt")


### filter min depth and alternative count
before_filter <- nrow(data)

min_depth <- 100
min_alt_count <- 3
data <- data[depth>min_depth & n_alt>min_alt_count]


after_filter <- nrow(data)
reduction <- round((before_filter-after_filter) / before_filter * 100, 2)

cat(paste(Sys.time(), "| min_depth and n_alt filtering |", 
          before_filter, "\t", after_filter, "\t", reduction, "\n"),
    file = log_file, append = TRUE)



### multiallelic removal
before_filter <- nrow(data)


# start is unqiue
data <- data |> 
  group_by(start) |> 
  filter(n() == 1) |> 
  ungroup()

# end is unqiue
data <- data |> 
  group_by(end) |> 
  filter(n() == 1) |> 
  ungroup()


after_filter <- nrow(data)
reduction <- round((before_filter-after_filter) / before_filter * 100, 2)

cat(paste(Sys.time(), "| Multi-allelic filtering |", 
          before_filter, "\t", after_filter, "\t", reduction, "\n"),
    file = log_file, append = TRUE)




### Somatic filter
#before_filter <- nrow(data)


#data <- data |> mutate(group = ifelse(n_alt < somatic_threshold, "somatic", "germline"))
#setDT(data)
#data <- data[group == "somatic"]


#after_filter <- nrow(data)
#reduction <- round((before_filter-after_filter) / before_filter * 100, 2)

#cat(paste(Sys.time(), "| Somatic/germline filtering |", 
#          before_filter, "\t", after_filter, "\t", reduction, "\n"),
#    file = log_file, append = TRUE)




### max depth: mean(depth) + 3*sd(depth)
# or 2*mean(depth)?
before_filter <- nrow(data)


max_depth <- mean(data$depth)+3*sd(data$depth)
setDT(data)
data <- data[depth<max_depth]


after_filter <- nrow(data)
reduction <- round((before_filter-after_filter) / before_filter * 100, 2)

cat(paste(Sys.time(), "| Max_depth filtering |", 
          before_filter, "\t", after_filter, "\t", reduction, "\n"),
    file = log_file, append = TRUE)




#write the file
write.table(data, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

