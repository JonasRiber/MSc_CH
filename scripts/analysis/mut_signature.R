

### author: Jonas Riber Jørgensen
### date: 22-04-2025
### desc: Combine the marked variants into one dataset
### output: mutational spectrum plots of all the variants together

library(tidyverse)
library(data.table)


setwd("/path/to/working/directory")
if(!interactive()) pdf(NULL) # suppress Rplots.pdf file

### setup args
args <- commandArgs(trailingOnly = TRUE)

#check if a file path was provided
if (length(args) < 2) {
  stop("Usage: Rscript mut_signature.R <input_directory.tsv> <output_directory>") # input: takes all files in the directory that ends with _full.tsv
}

folder_path <- args[1]   #First argument is the file path
#results/annotated_samples/
output_dir <- args[2]  #Second argument is output directory path
#results/all_variants_mutsig/

#check if in and out location exists
if (!file.exists(folder_path)) {
  stop(paste("File not found:", folder_path))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



#local test data
#folder_path <- "C:\\Users\\au683396\\Documents\\"





#list of all .tsv files
tsv_files <- list.files(folder_path, pattern = "\\_full.tsv$", full.names = TRUE)

nr_files <- length(tsv_files)



#function to read each file and add a column with the filename
read_and_tag <- function(file) {
  read.delim(file, stringsAsFactors = FALSE) %>%
    mutate(source_file = basename(file))
}

#apply the function to all files and combine them
combined_data <- bind_rows(lapply(tsv_files, read_and_tag))


combined_data <- combined_data |> filter(Alt_count<Somatic_threshold)


combined_data <- combined_data |> mutate(pt_id = sub("_.*", "", basename(source_file)))
id_translation <- read.table("/path/to/patientID_cipher", 
                             header = TRUE, sep = "\t", 
                             stringsAsFactors = F)
#local
#id_translation <- read.table("C:\\Users\\au683396\\Documents\\PatientIDs_paper.txt", header = TRUE, sep = "\t", stringsAsFactors = F)


combined_data <- merge(combined_data, id_translation, by.x = "pt_id", by.y = "ClusterID", all.x = T)


### overview table
overview_table <- combined_data |> group_by(status) |> 
  summarise(n = n(), 
            mean_VAF=mean(VAF), 
            mean_RD=mean(Total_count)) 
overview_table

output_table <- file.path(output_dir, "overview_table_variant_annotation.tsv")
write.table(overview_table, file = output_table, sep = "\t", quote = FALSE, row.names = FALSE)

summary(combined_data$VAF)
summary(combined_data$Total_count)

overview_table <- combined_data |> group_by(source_file, status) |> 
  summarise(n = n(), 
            mean_VAF=mean(VAF), 
            mean_RD=mean(Total_count)) 
overview_table

output_table <- file.path(output_dir, "overview_table_variant_annotation_per_patient.tsv")
write.table(overview_table, file = output_table, sep = "\t", quote = FALSE, row.names = FALSE)

# Vaf placement?
combined_data |> ggplot(aes(x = VAF, group = status, fill = status)) +
  geom_histogram(alpha=0.6) +
  NULL

# difference in depth?
combined_data |> ggplot(aes(x=Total_count, group = status)) +
  geom_histogram() +
  facet_grid(~status) +
  NULL



# large disparity between patients?
combined_data |> group_by(source_file, status) |> summarise(n=n())


# ratio between CH+ and not
status_ratio <- combined_data |> group_by(source_file, status, PaperID) |> 
  summarise(n=n(), .groups = "drop") |> 
  group_by(source_file) |> 
  mutate(ratio=n/sum(n),
         patient_id = sub("_.*", "", basename(source_file)))


# ratio of CH+ variants within each patient
status_ratio |> ggplot(aes(x = source_file, y = ratio, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Ratio by Source and Status",
       x = "Patient", y = "Ratio") +
  theme_minimal() +
  NULL

setDT(status_ratio)
true_ratio <- status_ratio |> 
  filter(status==TRUE) |> 
  mutate(pt_id = paste0("pt", row_number()))

# ratio comparison of CH+ variants within each patient
plot <- true_ratio |> 
  ggplot(aes(x = reorder(PaperID, -ratio), y = ratio, fill = status)) +
  geom_bar(stat = "identity") +
  labs(title = "Ratio for CH+",
       x = "Patient", y = "Ratio") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1)) +
  NULL
plot

plot_name <- file.path(output_dir, "all_variants_ratio_CHpos.png")
ggsave(plot_name, plot, width = 8, height = 4)

### Variants found within genes
df <- combined_data |> 
  filter(status == TRUE) 

df <- df |> mutate(patient_id = word(source_file, 1, sep = "_"))

variant_counter <- df |> group_by(patient_id) |> summarise(n=n())
print("Mean CH varaint count")
mean(variant_counter$n)

plot <- df |> ggplot(aes(x = fct_infreq(Gene))) +
  geom_histogram(fill = "cornflowerblue", stat = "count") +
  labs(x = "Gene",
       y = "n",
       caption = paste0("Variants found in the driverlist",
                        "\nNumber of patients: ", nr_files,
                        "\nNumber of Variants: ", nrow(df),
                        "\nNo patient grouping")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   hjust = 1)) +
  NULL
plot

plot_name <- file.path(output_dir, "Genes_in_matched_variants.png")
ggsave(plot_name, plot, width = 8, height = 4)

df |> group_by(patient_id) |> summarise(n=n()) 


### mutational signature


#normalize mutation direction to pyrimidine (C or T)
revcomp <- function(seq) {sapply(seq, function(s) {
  s_rev <- rev(utf8ToInt(s)) %>% 
    intToUtf8()  # reverse the string     
  chartr("ACGTacgt", "TGCAtgca", s_rev)      # complement   
}, USE.NAMES = FALSE) 
}

#create col with the substitution type
combined_data <- combined_data %>%
  mutate(
    status_label = ifelse(status, "CH+", "CH-"),
    substype = paste0(Ref, ">", Alt),
    norm_substype = ifelse(Ref %in% c("C", "T"),
                           paste0(Ref, ">", Alt),
                           paste0(revcomp(Ref), ">", revcomp(Alt)))
  )





mutation_freq <- combined_data |> 
  mutate(substitution = paste0(Ref, ">", Alt), 
         status_label = ifelse(status, "CH+", "CH-")) |> 
  filter(Ref %in% c("A", "C", "G", "T") & Alt %in% c("A", "C", "G", "T")) |> 
  group_by(status_label, norm_substype) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(status_label) |> 
  mutate(proportion = count / sum(count))

mutation_freq_counts <- mutation_freq |> group_by(status_label) |> summarise("CH+sum" = sum(count))
CH_pos_count <- mutation_freq_counts[[2]][1]
CH_neg_count <- mutation_freq_counts[[2]][2]

mutation_freq <- mutation_freq |> 
  mutate(status_label = factor(status_label, levels = c("CH+", "CH-")))

#plot
plot <- mutation_freq |> 
  ggplot(aes(x = norm_substype, y = proportion, fill = status_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CH-" = "#1f77b4", "CH+" = "#ff7f0e")) +
  labs(title = "Mutational Spectrum", x = "Substitution", y = "Proportion",
       fill = "CH Status",
       caption = paste0("Number of patients:", nr_files,
                        "\nCH-: ", CH_neg_count,"   CH+: ", CH_pos_count)) +
  theme_minimal() +
  NULL
plot

plot_name <- file.path(output_dir, "all_variants_mutational_spectrum.png")
ggsave(plot_name, plot, width = 5, height = 3)











### Trinucleotide context
# requires the surrounding context
# totaling 96 signatures
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("BSgenome", "BSgenome.Hsapiens.UCSC.hg38"))

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

BSgenome.Hsapiens.UCSC.hg38::Hsapiens[["chr1"]][50000:50010]
get_seq <- function(chrom, start, end) as.character(BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chrom]][start:end])
get_seq("chr1",500000, 500010)
# "AAGGTATCCTC"


combined_data <- combined_data |> 
  mutate(context = NA, 
         substype = NA,
         status_label = NA
         )


for (i in 1:nrow(combined_data)) {
  chrom <- combined_data$Chr[i]
  start <- combined_data$Start[i]
  end <- combined_data$End[i]
  
  buffer <- 1                     # How much to each side of the positions should be added
  ref_sequence <- get_seq(chrom, start-buffer, end+buffer)

  combined_data$context[i] <- ref_sequence
}



#normalize mutation direction to pyrimidine (C or T)
revcomp <- function(seq) {sapply(seq, function(s) {
  s_rev <- rev(utf8ToInt(s)) %>% 
    intToUtf8()  # reverse the string     
  chartr("ACGTacgt", "TGCAtgca", s_rev)      # complement   
  }, USE.NAMES = FALSE) 
}

combined_data <- combined_data %>%
  mutate(
    status_label = ifelse(status, "CH+", "CH-"),
    substype = paste0(Ref, ">", Alt),
    norm_substype = ifelse(Ref %in% c("C", "T"),
                      paste0(Ref, ">", Alt),
                      paste0(revcomp(Ref), ">", revcomp(Alt))),
    norm_context = ifelse(Ref %in% c("C", "T"),
                          context,
                          revcomp(context))
  )


#plot
trinucleotide_freq <- combined_data |> 
  filter(Ref %in% c("A", "C", "G", "T") & Alt %in% c("A", "C", "G", "T")) |> 
  group_by(status_label, norm_substype, norm_context) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(status_label) |> 
  mutate(proportion = count / sum(count))


substype_fill_var <- c("C>A"="turquoise", "C>G"="black", "C>T"="red",
                       "T>A"="grey", "T>C"="palegreen", "T>G"="pink")
plot <- trinucleotide_freq |> 
  ggplot(aes(norm_context, proportion, fill=norm_substype)) +
  geom_col() +
  scale_fill_manual(values = substype_fill_var, guide="none") +
  facet_grid(.~norm_substype, scales="free_x", space="free") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=.5, hjust=0)) +
  facet_grid(status_label~norm_substype, scales="free_x", space="free") + 
  labs(title = "Trinucleotide spectrum",
       x = "Context", y = "Proportion",
       caption = paste0("Number of patients:", nr_files,
                        "\nCounts CH-: ", CH_neg_count,"   CH+: ", CH_pos_count)) +
  NULL
plot

plot_name <- file.path(output_dir, "all_variants_trinucleotid_spectrums.png")
ggsave(plot_name, plot, width = 10, height = 7, dpi = 600)

# cosine similarity
#sort for mutation type
setDT(trinucleotide_freq)

cossim <- function(A,B) { (sum(A*B))/sqrt((sum(A^2))*(sum(B^2))) }
cosin_sim <- cossim(trinucleotide_freq[status_label=="CH+", proportion], trinucleotide_freq[status_label=="CH-", proportion])
cosin_sim <- format(cosin_sim, big.mark=",", digits=2)



### differences
tri_differences <- trinucleotide_freq |> select(status_label, norm_substype, norm_context, proportion) |> 
  pivot_wider(names_from = status_label,
              values_from = proportion,
              names_prefix = "prop_") |> 
  mutate(proportion_diff = `prop_CH+` - `prop_CH-`)

all_types_present <- any(is.na(tri_differences$`prop_CH+`) & is.na(tri_differences$`prop_CH-`))


plot <- tri_differences |> 
  ggplot(aes(norm_context, proportion_diff, fill=norm_substype)) +
  geom_col() +
  scale_fill_manual(values = substype_fill_var, guide="none") +
  facet_grid(.~norm_substype, scales="free_x", space="free") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=.5, hjust=0),
        plot.margin = margin(10,60,10,10),
        clip = "off") +
  labs(title = "Trinucleotide spectrum",
       x = "Context", y = "Proportion difference (CHpos - CHneg)",
       caption = paste0("Number of patients:", nr_files,
                        "\nCounts CH-: ", CH_neg_count,"   CH+: ", CH_pos_count,
                        "\nCosine simularity: ", cosin_sim,
                        "\nAll mut.types present: ", all_types_present)) + # Cosine simularity would require all substypes to be present
  NULL

plot_name <- file.path(output_dir, "all_variants_trinucleotide_diff_w_cossim.png")
ggsave(plot_name, plot, width = 10, height = 4, dpi = 600)







### potentially make spectrum for each patient?











