

### author: Jonas Riber Jørgensen
### date: 30-04-2025
### desc: plots of the fisher method tests and add FDR adjustment
### output: tsv file with patient and their CH status


setwd("/path/to/working/directory")

library(tidyverse)
library(data.table)

if(!interactive()) pdf(NULL) # suppress Rplots.pdf file

data_path <- "results/variant_table_test_results_var_ER.tsv"

#local test
#data_path <- "C:\\Users\\au683396\\Documents\\variant_table_test_results.tsv"


data <- read.delim(data_path, header = TRUE, sep = "\t")
nr_patients <- nrow(data) 


threshold <- 0.05 # pvalue threshold
nr_CH_patients <- sum(data$combined_pvalue<0.05)

plot <- data |> ggplot(aes(x = combined_pvalue)) +
  geom_histogram(aes(fill = ..x.. < threshold),
                 bins = 30) +
  scale_fill_manual(name = NULL,
                    values = c("TRUE" = "cornflowerblue", "FALSE" = "firebrick"),
                    labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold")) +
  labs(title = "Distribution of fused pvalues",
       x = "Fused pvalues",
       y = "n",
       caption = paste0("Fused pvalue distribution, the threshold is set to 0.05",
                        "\nNr of patients: ", nr_patients,
                        "\nNr of sig: ", nr_CH_patients)) +
  theme_minimal() +
  NULL
plot

plot_name <- file.path("results/plots/pvalue_distribution.png")
ggsave(plot_name, plot, width = 8, height = 3, units = "in", dpi = 300)



data |> filter(combined_pvalue < threshold) |> arrange(combined_pvalue)




### Use FDR
data$adjusted_pval <- p.adjust(data$combined_pvalue, method = "BH")
# done to control the expected proportion of false positives amoing the results we call sig
# with FDR at 0.05 it means that on average no more than 5% of the hypothesis
# we call sig are expected to be false positives 

BH_threshold <- 0.1
CH_patients <- data |> filter(adjusted_pval<BH_threshold) |>  arrange(adjusted_pval)
CH_patients
nr_CH_patients <- nrow(CH_patients)

percentage_CH <- nr_CH_patients/nr_patients * 100
nr_CH_patients
percentage_CH


plot <- data |> ggplot(aes(x = adjusted_pval)) +
  geom_histogram(aes(fill = ..x.. < BH_threshold),
                 bins = 50) + 
  scale_fill_manual(name = NULL,
                    values = c("TRUE" = "cornflowerblue", "FALSE" = "firebrick"),
                    labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold")) +
  labs(title = "pvalue distribution",
       x = "pvalue (FDR adjusted)",
       y = "n", 
       caption = paste0("Number of patients: ", nrow(data),
                        "\nNumber of sig: ", nr_CH_patients,
                        "\nBH threshold: ", BH_threshold)) +
  theme_minimal() +
  NULL
plot

plot_name <- file.path("results/plots/FDR_adjusted_distribution.png")
ggsave(plot_name, plot, width = 8, height = 3, units = "in", dpi = 300)







# write file with the patient status
data <- data |> mutate(patient_status = ifelse(data$patient_id %in% CH_patients$patient_id, "positive", "negative"))


write.table(x = data, file = "results/patient_status_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE)





#df_path <- "C:\\Users\\au683396\\Documents\\variant_table_mms_depth.tsv"


#df <- read.delim(df_path, header = TRUE, sep = "\t")




# --- redo mutational signature --- #
# with new grouping
#folder_path <- "C:\\Users\\au683396\\Documents\\"
output_dir <- "results/patient_mutsig/"
folder_path <- "results/annotated_samples/"
#list of all .tsv files

tsv_files <- list.files(folder_path, pattern = "\\_full.tsv$", full.names = TRUE)

nr_files <- length(tsv_files)



#function to read each file and add a column with the filename
read_and_tag <- function(file) {
  id <- sub("_.*", "", basename(file))
  read.delim(file, stringsAsFactors = FALSE) %>%
    mutate(patient_id = id
           )
}

#apply the function to all files and combine them
combined_data <- bind_rows(lapply(tsv_files, read_and_tag))
combined_data <- combined_data |> filter(Alt_count<Somatic_threshold)
combined_data <- combined_data |> filter(Total_count>100 & Alt_count>3)

combined_data <- combined_data |> 
  mutate(patient_CH_status = ifelse(patient_id %in% CH_patients$patient_id, "positive", "negative"))



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


mutation_freq <- combined_data |> 
  mutate(substitution = paste0(Ref, ">", Alt), 
         status_label = ifelse(status, "CH+", "CH-"),
         group_label = ifelse(patient_CH_status == "positive" & status_label == "CH+",
                              "CH+ pt/drivers",
                              "Other")) |> 
  filter(Ref %in% c("A", "C", "G", "T") & Alt %in% c("A", "C", "G", "T")) |> 
  group_by(group_label, norm_substype) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(group_label) |> 
  mutate(proportion = count / sum(count))

mutation_freq_counts <- mutation_freq |> group_by(group_label) |> summarise("CH+sum" = sum(count))
CH_pos_count <- mutation_freq_counts[[2]][1]
CH_neg_count <- mutation_freq_counts[[2]][2]
CH_pos_count
CH_neg_count

#plot
plot <- mutation_freq |> 
  ggplot(aes(x = norm_substype, y = proportion, fill = group_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Other" = "#1f77b4", "CH+ pt/drivers" = "#ff7f0e")) +
  labs(title = "Mutational Spectrum", x = "Substitution", y = "Proportion",
       fill = "CH Status",
       caption = paste0("Number of patients:", nr_files,
                        "\nCH-: ", CH_neg_count,"   CH+: ", CH_pos_count)) +
  theme_minimal() +
  NULL
plot


plot_name <- file.path(output_dir, "significant_patients_mutational_signature.png")
ggsave(plot_name, plot, width = 5, height = 3, dpi = 600)





### tri-nucleotide ###

trinucleotide_freq <- combined_data |> 
  filter(Ref %in% c("A", "C", "G", "T") & Alt %in% c("A", "C", "G", "T")) |> 
  group_by(patient_CH_status, norm_substype, norm_context) |> 
  summarise(count = n(), .groups = "drop") |> 
  group_by(patient_CH_status) |> 
  mutate(proportion = count / sum(count))


CH_neg_count <- trinucleotide_freq |> 
  filter(patient_CH_status=="negative") |> 
  summarise(n = sum(count))
CH_neg_count <- CH_neg_count[[2]][1]


CH_pos_count <- trinucleotide_freq |> 
  filter(patient_CH_status=="positive") |> 
  summarise(n = sum(count))
CH_pos_count <- CH_pos_count[[2]][1]


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
  facet_grid(patient_CH_status~norm_substype, scales="free_x", space="free") + 
  labs(title = "Trinucleotide spectrum",
       x = "Context", y = "Proportion",
       caption = paste0("Number of patients:", nr_files,
                        "\nCH- patients: ", nr_patients-nr_CH_patients, "   CH+ patients: ", nr_CH_patients,
                        "\nMms in neg group: ", CH_neg_count,"   Mms in pos group: ", CH_pos_count)) +
  NULL
plot

plot_name <- file.path(output_dir, "patient_level_trinucleotid_spectrums.png")
ggsave(plot_name, plot, width = 10, height = 7, dpi = 600)








# cosine similarity
#sort for mutation type
setDT(trinucleotide_freq)

cossim <- function(A,B) { (sum(A*B))/sqrt((sum(A^2))*(sum(B^2))) }
cosin_sim <- cossim(trinucleotide_freq[patient_CH_status=="positive", proportion], trinucleotide_freq[patient_CH_status=="negative", proportion])
cosin_sim <- format(cosin_sim, big.mark=",", digits=2)



### differences
tri_differences <- trinucleotide_freq |>
  select(patient_CH_status, norm_substype, norm_context, proportion) |>
  pivot_wider(names_from = patient_CH_status, 
              values_from = proportion) |>
  mutate(proportion_diff = positive - negative)

all_types_present <- any(is.na(tri_differences$positive) & is.na(tri_differences$negative))


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
       caption = paste0("Number of patients:", nr_patients,
                        "\nCounts CH-: ", CH_neg_count,"   CH+: ", CH_pos_count,
                        "\nCosine simularity: ", cosin_sim,
                        "\nAll mut.types present: ", all_types_present)) + # Cosine simularity would require all substypes to be present
  NULL

plot_name <- file.path(output_dir, "patient_level_trinucleotide_diff_w_cossim.png")
ggsave(plot_name, plot, width = 10, height = 7, dpi = 600)



