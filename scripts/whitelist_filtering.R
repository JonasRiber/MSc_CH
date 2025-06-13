

### author: Jonas Riber Jørgensen
### date: 20-04-2025
### desc: based on the whitelist collected from the blood article (which includes the variants from IntOGen) 
#         mark the status if a variant is found within the whitelist. 
### output: 2 files, first for the full input of variants marked with T/F in new col called status
#           Second file is filtered to only include the Ch variants 






### Initialize and load data

library(tidyverse)
library(data.table)
setwd("/Path/to/working/directory")

if(!interactive()) pdf(NULL) # suppress Rplots.pdf file


### setup args
args <- commandArgs(trailingOnly = TRUE)

#check if a file path was provided
if (length(args) < 2) {
  stop("Usage: Rscript whitelist_filtering.R <input_file.tsv> <output_directory>")
}

file_path <- args[1]   #First argument is the file path
output_dir <- args[2]  #Second argument is output directory path


#check if in and out location exists
if (!file.exists(file_path)) {
  stop(paste("File not found:", file_path))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


#read input
data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(data)
colnames(data)

#get id
file_name <- tools::file_path_sans_ext(basename(file_path))
id <- print(str_extract(string = file_name, pattern = "PT_ID_PATTERN"))

#load whitelists
gList<-fread("backup/whitelists/Full_CHIP_gene_list_08262022.txt")
whitelist.mis<-fread("backup/whitelists/CHIP_missense_vars_cv_08262022.txt")
whitelist.splice<-fread("backup/whitelists/CHIP_splice_vars_cv_08262022.txt")
whitelist.nons_and_FS<-fread("backup/whitelists/CHIP_nonsense_FS_vars_cv_08262022.txt")


setDT(data)




### setup CH status
data <- data |> mutate(status = F)

#remove variants above somatic thrshold (was kept for prior plotting)
data <- data |> filter(Alt_count<Somatic_threshold)


### Specific variant seen in driverlist - based on AAchange
data$Gene <- as.character(data$Gene)
data$AA_change <- as.character(data$AA_change)
whitelist.mis$Gene <- as.character(whitelist.mis$Gene)
whitelist.mis$AAChange <- as.character(whitelist.mis$AAChange)


for (i in 1:nrow(data)) {
  gene <- data$Gene[i]
  variant_string <- data$AA_change[i]
  aa_entries <- unlist(strsplit(variant_string, ","))
  aa_changes <- sub(".*:p\\.([A-Z0-9]+)$", "\\1", aa_entries)

  match <- F
  for (aa_change in aa_changes) {

    #check pair of gene and any AAchange matches wgitelist
    if (any(whitelist.mis$Gene == gene & whitelist.mis$AAChange == aa_change)){
      match <- T
      print("Match!")
      break
    }
  }
  
  data$status[i] <- match
}


### nonsyn
for (i in 1:nrow(data)){
  if (data$Gene[i] %in% gList$Gene && data$Exonic_func[i] == "nonsynonymous SNV") {
    data$status[i] <- TRUE
  }
} 


### stopgain
for (i in 1:nrow(data)){
  if (data$Gene[i] %in% whitelist.nons_and_FS$Gene && data$Exonic_func[i] == "stopgain") {
    data$status[i] <- TRUE
  }
} 


### splice site variants
for (i in 1:nrow(data)){
  if (data$Gene[i] %in% whitelist.splice$Gene && data$Func[i] == "splicing") {
    data$status[i] <- TRUE
  }
} 


### non-frameshift multinuecleotide substitutions

#SF3A1:NM_005877.6:exon13:c.G2010T:p.V670V
#TEX13D:NM_001355534.2:exon1:c.1086_1087delinsCC:p.M362_D363delinsIH


### exceptions
# ASXL1 and ASXL2 - Frameshift/nonsense/splice-site in exon 11-12
exception_list <- c("ASXL1", "ASXL2")
exception_critiria <- c("exon11", "exon12")

for (i in 1:nrow(data)) {
  if (data$Gene[i] %in% exception_list) {
    variant_string <- data$AA_change[i]
    exons <- str_extract_all(variant_string, "exon\\d+")[[1]]
    
    for (exon in exons) {
      if (exon %in% exception_critiria && 
          data$Exonic_func[i]=="splicing" || 
          data$Exonic_func[i]=="nonsynonymous SNV" || 
          data$Exonic_func[i]=="frameshift") {
        data$status[i] <- T
        break
      }
    }
  }
}



tmp <- data[Gene=="CBLB" & Func=="exonic"]

# CBL - RING finger missense p.381-421
exception_criteria <- c(381, 421)

for (i in 1:nrow(data)) {
  if (data$Gene[i] == "CBL" && data$Func[i] == "exonic") {
    variant_string <- data$AA_change[i]
    protein_change <- str_extract_all(variant_string, "p\\.[^,]+")[[1]]
    protein_positions <- str_extract(protein_change, "\\d+")
    
    for (pos in protein_positions) {
      if (pos>exception_criteria[1] &&
          pos<exception_criteria[2] &&
          data$Exonic_func[i]=="nonsynonymous SNV") {
        data$status[i] <- T
        break
      }
    }
  }
}



# CBLB - RING finger missense p.372-412
exception_criteria <- c(372, 412)

for (i in 1:nrow(data)) {
  if (data$Gene[i] == "CBLB" && data$Func[i] == "exonic") {
    variant_string <- data$AA_change[i]
    protein_change <- str_extract_all(variant_string, "p\\.[^,]+")[[1]]
    protein_positions <- str_extract(protein_change, "\\d+")

    print(protein_positions)
    for (pos in protein_positions) {
      if (pos>exception_criteria[1] &&
          pos<exception_criteria[2] &&
          data$Exonic_func[i]=="nonsynonymous SNV") {
          
        data$status[i] <- T
          break
        }
     }
  }
}


# PPM1D Frameshift/nonsense,  exon 5 or 6
exception_list <- c("PPM1D")
exception_critiria <- c("exon5", "exon6")

for (i in 1:nrow(data)) {
  if (data$Gene[i] %in% exception_list) {
    variant_string <- data$AA_change[i]
    exons <- str_extract_all(variant_string, "exon\\d+")[[1]]
    
    for (exon in exons) {
      if (exon %in% exception_critiria && 
          data$Exonic_func[i]=="splicing" || 
          data$Exonic_func[i]=="nonsynonymous SNV" || 
          data$Exonic_func[i]=="frameshift") {
        data$status[i] <- T
        break
      }
    }
  }
}


# ZBTB33 - missense or in-frame indel in catalytic domains (p.9-126, 332-591)
exception_criteria1 <- c(9, 126)
exception_criteria2 <- c(332, 591)

for (i in 1:nrow(data)) {
  if (data$Gene[i] == "ZBTB33" && data$Func[i] == "exonic") {
    variant_string <- data$AA_change[i]
    protein_change <- str_extract_all(variant_string, "p\\.[^,]+")[[1]]
    protein_positions <- str_extract(protein_change, "\\d+")
    
    print(protein_positions)
    for (pos in protein_positions) {
      if (pos>exception_criteria1[1] &&
          pos<exception_criteria1[2] &&
          data$Exonic_func[i]=="nonsynonymous SNV") {
        
        data$status[i] <- T
        break
      }
      
      if (pos>exception_criteria2[1] &&
          pos<exception_criteria2[2] &&
          data$Exonic_func[i]=="nonsynonymous SNV") {
        
        data$status[i] <- T
        break
      }
    }
  }
}





# TET2 -  missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
#handled by nonsyn mut


# CSF3R - truncating c.741-791 
# handled by stopgain





### write marked file
name <- file.path(output_dir, paste0(id, "_CH_marked_full.tsv"))
write.table(data, name, sep = "\t", quote = FALSE, row.names = FALSE)

### write filtered version
filtered_data <- data |> filter(status==T)    

name <- file.path(output_dir, paste0(id, "_CH_marked_filtered.tsv"))
write.table(filtered_data, name, sep = "\t", quote = FALSE, row.names = FALSE)








