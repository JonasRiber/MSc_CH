




### author: Jonas Riber Jørgensen
### date: 24-04-2025
### desc: combine the CH+ annotated variants into a table with var in rows and patient in col
### output: outputs summary tables with the mms | depth for each patient per variant and the results of fishers method per patient

library(tidyverse)
library(data.table)
if(!interactive()) pdf(NULL) # suppress Rplots.pdf file

setwd("/Path/to/working/directory")


variant_dir <- "results/annotated_samples/" 
zero_dir <- "data/0_positions/"
depth_dir <- "data/0_positions/merged_depths/"



# Read all depth files
depth_files <- list.files(depth_dir, pattern = "merged_0pos_depth.tsv$", full.names = TRUE)

depth_data <- lapply(depth_files, function(file) {
  patient_id <- str_extract(basename(file), "^[^_]+")
  read_tsv(file, col_names = c("Chr", "Start", "End", "Ref", "Depth"), col_types = cols()) |> 
    mutate(patient_id = patient_id)
}) |> bind_rows()



# Read the 0-pos files
zero_files <- list.files(zero_dir, pattern = "_0_positions.tsv$", full.names = TRUE)

zero_data <- lapply(zero_files, function(file) {
  patient_id <- str_extract(basename(file), "^[^_]+")
  read_tsv(file, col_names = TRUE, col_types = cols()) |> 
    rename(patient_id = id) |> 
    mutate(patient_id = patient_id)
}) |> bind_rows()



# match depths to positions
zero_data <- left_join(zero_data, depth_data, by = c("Chr", "Start", "End", "Ref", "patient_id"))

#sum(is.na(zero_data$Depth))
# all entries should be filled so no NAs


#adjust for the indexing change for depth collection
zero_data$Start <- zero_data$Start+1


#reformat splicing sites
zero_data$AA_change[zero_data$Func == "splicing"] <- paste0(
  zero_data$Gene[zero_data$Func=="splicing"], ":",
  zero_data$Start[zero_data$Func== "splicing"], "-",
  zero_data$End[zero_data$Func== "splicing"], ":splicing"
)




# Re-gather all data files
tsv_files <- list.files(variant_dir, pattern = "\\marked_filtered.tsv$", full.names = TRUE)


# function to read each file and add a column with the id from filename
#read_and_tag <- function(file) {
#  read_tsv(file, col_types = cols(.default = "c")) |> 
#    mutate(id = strsplit(basename(file), "_")[[1]][1])
#}


read_and_tag <- function(file) {
  df <- read_tsv(file,
                 col_types = cols(
                   Chr = col_character(),
                   Start = col_double(),
                   End = col_double(),
                   Ref = col_character(),
                   Alt = col_character(),
                   Func = col_character(),
                   Gene = col_character(),
                   Exonic_func = col_character(),
                   AA_change = col_character(),
                   Gnomad_all = col_character(),
                   Clinvar = col_character(),
                   Cosmic = col_character(),
                   Mtype = col_double(),
                   Cons_count = col_double(),
                   Alt_count = col_double(),
                   Total_count = col_double(),
                   VAF = col_double(),
                   Somatic_threshold = col_double(),
                   status = col_character()
                 ))
  
  df <- df |> mutate(id = strsplit(basename(file), "_")[[1]][1])
  
  return(df)
}

# apply the function to all files and combine them
data <- bind_rows(lapply(tsv_files, read_and_tag))
data <- data |> filter(Alt_count<Somatic_threshold)
data <- data |> rename(patient_id = id, Depth = Total_count)

data$AA_change[data$Func == "splicing"] <- paste0(
  data$Gene[data$Func=="splicing"], ":",
  data$Start[data$Func== "splicing"], "-",
  data$End[data$Func== "splicing"], ":splicing"
)





# combine variant and zero_pos
df <- bind_rows(data, zero_data) |> 
  distinct(Chr, Start, End, Ref, AA_change, patient_id, .keep_all = TRUE)
df$Alt_count[is.na(df$Alt_count)] <- 0





### re-organize table
df <- df |> mutate(mm_depth = paste0(Alt_count, " | ", Depth))

df_wide <- df |> 
  select(AA_change, patient_id, mm_depth) |> 
  pivot_wider(
    names_from = patient_id,
    values_from = mm_depth
  )





### remove odd variant - SETBP1
#df_wide <- df_wide |> filter(AA_change!="SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T")



# sum mms in each row 
# sum depth in each row
# divide these

df_wide <- df_wide |> 
  rowwise() |> 
  mutate(
    #separate numbers in all columns except the first one
    sum_mms = sum(as.numeric(sapply(c_across(2:ncol(df_wide)), function(x) strsplit(x, " | ")[[1]][1]))),
    sum_depths = sum(as.numeric(sapply(c_across(2:ncol(df_wide)), function(x) strsplit(x, " | ")[[1]][3]))),
    error_rate = sum_mms / sum_depths
  ) %>%
  ungroup()


# find mean of error rate
p_hat <- mean(df_wide$error_rate) 
print("The mean error rate across all variants:")
print(p_hat)



### save the current table with mms | depth in each cell
write.table(df_wide, "results/variant_table_mms_depth.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



#--- mean ER ---#
# binomial test in each cell with the p_hat for mean ER
binomial_test_from_string <- function(cell) {
  #extract mismatches and depth from the "mismatches | depth"
  mms <- as.numeric(strsplit(cell, " | ")[[1]][1])
  depth <- as.numeric(strsplit(cell, " | ")[[1]][3])
  
  #run binomial test 
  result <- binom.test(mms, depth, p = p_hat)
  return(result$p.value)
}


# patient cols
patient_columns <- 2:(ncol(df_wide) - 3) 

#apply the function to only the patient columns and replace the mismatch-depth strings with p-values
df_result <- df_wide
df_result[, patient_columns] <- apply(df_result[, patient_columns], c(1, 2), function(cell) binomial_test_from_string(cell))




#Fishers method
#init dataframe for results
results_df <- data.frame(
  patient_id = character(),
  chisq_statistic = numeric(),
  combined_pvalue = numeric(),
  stringsAsFactors = FALSE
)
setDT(results_df)

for (i in patient_columns) {
  current_patient <- df_result[i]
  patient_id <- colnames(df_result[i])
  
  # calc combine test stat
  chi_sq <- -2*sum(log(current_patient))
  
  #deg of freedom (2x pvalues)
  deg_free <- 2 * nrow(current_patient)
  
  #use X2 stat to get combined pval
  combi_pval <- pchisq(q = chi_sq, df = deg_free, lower.tail = FALSE)
  
  results_df <- rbind(results_df, data.frame(patient_id = patient_id,
                               chisq_statistic = chi_sq,
                               combined_pvalue = combi_pval,
                               stringsAsFactors = FALSE)
        )
}

print("degrees of freedom used")
print(deg_free)


write.table(results_df, "results/variant_table_test_results_mean_ER.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



results_df$adjusted_pval <- p.adjust(results_df$combined_pvalue, method = "BH")

BH_threshold <- 0.1

obs_sig <- sum(results_df$adjusted_pval<BH_threshold)
print("Number of significant patients with mean ER")
obs_sig













#--- Varient ER ---#
# binomial test in each cell with the ER per variant 
binomial_test_ER_per_var <- function(cell, er) {
  #extract mismatches and depth from the "mismatches | depth"
  mms <- as.numeric(strsplit(cell, " | ")[[1]][1])
  depth <- as.numeric(strsplit(cell, " | ")[[1]][3])

  
  #run binomial test 
  result <- binom.test(mms, depth, p = er)
  #print(result$p.value)
  return(result$p.value)
}


# patient cols
patient_columns <- 2:(ncol(df_wide) - 3) 

#apply the function to only the patient columns and replace the mismatch-depth strings with p-values
df_var_ER <- df_wide

for (i in 1:nrow(df_var_ER)){
  er <- as.numeric(df_var_ER$error_rate[i])
  for (j in patient_columns){
    current_cell <- df_var_ER[[j]][i]
    pval <- binomial_test_ER_per_var(current_cell, er)
    
    df_var_ER[[j]][i] <- pval
    
  }
} 

df_var_ER[df_var_ER == TRUE] <- NA


write.table(df_var_ER, "results/variant_table_var_ER_pvalues.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#df_var_ER <- na.omit(df_var_ER)

df_var_ER[, patient_columns] <- lapply(df_var_ER[, patient_columns], function(x) as.numeric(as.character(x)))


#Fishers method
#init dataframe for results
results_two_df <- data.frame(
  patient_id = character(),
  chisq_statistic = numeric(),
  combined_pvalue = numeric(),
  stringsAsFactors = FALSE
)
setDT(results_two_df)

for (i in patient_columns) {
  current_patient <- df_var_ER[i]
  patient_id <- colnames(df_var_ER[i])

  # calc combine test stat
  chi_sq <- -2*sum(log(current_patient))
    
  #deg of freedom (2x pvalues)
  deg_free <- 2 * nrow(current_patient)
  
  #use X2 stat to get combined pval
  combi_pval <- pchisq(q = chi_sq, df = deg_free, lower.tail = FALSE)
  
  results_two_df <- rbind(results_two_df, data.frame(patient_id = patient_id,
                                             chisq_statistic = chi_sq,
                                             combined_pvalue = combi_pval,
                                             stringsAsFactors = FALSE)
  )
}

print("degrees of freedom used")
print(deg_free)


write.table(results_two_df, "results/variant_table_test_results_var_ER.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



results_two_df$adjusted_pval <- p.adjust(results_two_df$combined_pvalue, method = "BH")

obs_sig <- sum(results_two_df$adjusted_pval<BH_threshold)
print("Number of significant patients with variant level ER")
obs_sig



#--- subsample - k-fold crossvalidation ---#

combine_pvals_fisher <- function(df, patient_columns) {
  #initialize output
  fisher_result <- data.frame(
    patient_id = character(),
    chisq_statistic = numeric(),
    combined_pvalue = numeric(),
    p_adjusted = numeric(),
    stringsAsFactors = FALSE
  )
  setDT(fisher_result)
  
  
  for (i in patient_columns) {
    current_patient <- df[i]
    patient_id <- colnames(df[i])
    
    # Fishers method
    chi_sq <- -2*sum(log(current_patient))
    deg_free <- 2 * nrow(current_patient)
    combi_pval <- pchisq(q = chi_sq, df = deg_free, lower.tail = FALSE)
    
    #FDR adjust
    p_adjusted <- p.adjust(p = combi_pval, method = "BH")
    
    #store result
    fisher_result <- rbind(fisher_result, data.frame(patient_id = patient_id,
                                                     chisq_statistic = chi_sq,
                                                     combined_pvalue = combi_pval,
                                                     p_adjusted = p_adjusted,
                                                     stringsAsFactors = FALSE)
    )
  }
  return(fisher_result)
}


rows <- sample(nrow(df_var_ER))
half <- floor(nrow(df_var_ER) / 2)

df1 <- df_var_ER[rows[1:half], ]
df2 <- df_var_ER[rows[(half + 1):nrow(df_var_ER)], ]

df1_fisher <- combine_pvals_fisher(df1, patient_columns)
df1_fisher$split <- "A"

df2_fisher <- combine_pvals_fisher(df2, patient_columns)
df2_fisher$split <- "B"

df_subsampling <- rbind(df1_fisher, df2_fisher)
df_subsampling[p_adjusted<BH_threshold]


sig_A <- df1_fisher[p_adjusted < BH_threshold, "patient_id"]
sig_B <- df2_fisher[p_adjusted < BH_threshold, "patient_id"]

common_sig_pt <- intersect(sig_A, sig_B)


print("With a 50:50 split, the following patients were significant in both")
print(common_sig_pt)
print(paste0("Totaling: ", length(common_sig_pt[[1]])))


write.table(df_subsampling, "results/subsampling_variant_table_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)




#--- Repeat resampling ---#

# Initialize storage
overlap_counts <- numeric(100)
n_sig_A <- numeric(100)
n_sig_B <- numeric(100)
patient_overlap_list <- list()

for (i in 1:100) {
  rows <- sample(nrow(df_var_ER))
  half <- floor(nrow(df_var_ER) / 2)
  
  df1 <- df_var_ER[rows[1:half], ]
  df2 <- df_var_ER[rows[(half + 1):nrow(df_var_ER)], ]
  
  df1_fisher <- combine_pvals_fisher(df1, patient_columns)
  df2_fisher <- combine_pvals_fisher(df2, patient_columns)
  
  sig_A <- df1_fisher[df1_fisher$p_adjusted < BH_threshold, "patient_id"]
  sig_B <- df2_fisher[df2_fisher$p_adjusted < BH_threshold, "patient_id"]
  
  common <- intersect(sig_A, sig_B)
  
  # Record stats
  n_sig_A[i] <- length(sig_A)
  n_sig_B[i] <- length(sig_B)
  overlap_counts[i] <- length(common)
  patient_overlap_list[[i]] <- common
}

# Summary stats
cat("Mean overlap:", mean(overlap_counts), "\n")
cat("Median overlap:", median(overlap_counts), "\n")
cat("Mean # of significant patients in split A:", mean(n_sig_A), "\n")
cat("Mean # of significant patients in split B:", mean(n_sig_B), "\n")

# Optional: frequency table of how often each patient appeared in the overlap
all_common_patients <- unlist(patient_overlap_list)
common_freq_table <- sort(table(all_common_patients), decreasing = TRUE)

# View top consistently selected patients
head(common_freq_table)

all_common_patients <- unlist(patient_overlap_list)
patient_freq <- as.data.frame(table(all_common_patients))
colnames(patient_freq) <- c("patient_id", "frequency")

patient_freq <- patient_freq[order(-patient_freq$frequency), ]


id_translation <- read.table("~/path/to/patient_id_cipher", 
                             header = TRUE, sep = "\t", 
                             stringsAsFactors = F)

patient_freq <- merge(patient_freq, id_translation, by.x = "patient_id", by.y = "ClusterID", all.x = T)


plot <- ggplot(patient_freq, aes(x = reorder(PaperID, -frequency), y = frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Frequency of Patient Selection Across Resampling",
       x = "Patient ID",
       y = "Times Selected in Both Splits") +
  NULL
plot


plot_name <- file.path("results/variant_table_plots/Resampling_patient_count.png")
ggsave(plot_name, plot, width = 6, height = 4)








#--- Permutation testing - shuffling labels ---#
# is the observed-per-patient p-values significantly different from what we can expect by random chance

n_perm <- 10000
permute_results <- data.frame(
  iteration = integer(),
  sig_patients = integer(),
  patient_ids = character(),
  stringsAsFactors = FALSE
)


### shuffle rows independently
shuffle_selected_cols <- function(df) {
  cols_to_shuffle <- 2:(ncol(df) - 3)
  df[cols_to_shuffle] <- t(apply(df[cols_to_shuffle], 1, sample))
  return(df)
}


### shuffle everything
#shuffle_selected_cols <- function(df) {
#  cols_to_shuffle <- 2:(ncol(df)-3)
#  values <- unlist(df[cols_to_shuffle])
#  shuffled_values <- sample(values)
#  
#  reshaped <- matrix(shuffled_values, nrow= nrow(df))
#  df[cols_to_shuffle] <- as.data.frame(reshaped)
#  
#  return(df)
#}



print("Permutation testing:")

for (i in 1:n_perm) {
  shuffled_sample <- shuffle_selected_cols(df_var_ER)
  
  
  #Fishers method
  #init dataframe for results
  x <- data.frame(
    patient_id = character(),
    chisq_statistic = numeric(),
    combined_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  setDT(x)
  
  for (j in patient_columns) {
    current_patient <- shuffled_sample[j]
    patient_id <- colnames(shuffled_sample[j])
    
    # calc combine test stat
    chi_sq <- -2*sum(log(current_patient))
    
    #deg of freedom (2x pvalues)
    deg_free <- 2 * nrow(current_patient)
    
    #use X2 stat to get combined pval
    combi_pval <- pchisq(q = chi_sq, df = deg_free, lower.tail = FALSE)
    
    x <- rbind(x, data.frame(patient_id = patient_id,
                             chisq_statistic = chi_sq,
                             combined_pvalue = combi_pval,
                             stringsAsFactors = FALSE)
    )
  }
  
  
  #FDR adjust
  x$adjusted_pval <- p.adjust(x$combined_pvalue, method = "BH")
  
  selected <- x[x$adjusted_pval < BH_threshold, ]
  selected_ids <- unique(selected$patient_id)
  
  #collect sig patentis
  permute_results <- rbind(permute_results, data.frame(
    iteration = i,
    sig_patients = length(selected_ids),
    patient_ids = paste(selected_ids, collapse = ";"))
  )
}


write.table(permute_results, "results/permutations_of_variant_table_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


permuted <- sum(permute_results$sig_patients>=obs_sig)
emp_p <- permuted/n_perm


print(paste0("Out of ", n_perm, " permutations, ",
             permuted, " were either as or more extreme"))
print(paste0("empirical p-value: ", emp_p))


permute_plot <- permute_results |> 
  ggplot(aes(x = sig_patients)) +
  geom_histogram(fill = "firebrick", binwidth = 1, 
                 boundary = 0.5, ) +
  geom_vline(xintercept = obs_sig, linetype = "dashed", 
             show.legend = TRUE, linewidth = 1,
             color = "cornflowerblue") +
  labs(title = "Null distribution of permutation test", 
       x = "Significant patients",
       y = "n",
       caption = paste0("Number of permutations: ", n_perm,
                        "\nobserved significant: ", obs_sig,
                        "\npermuted as or more extreme: ", permuted)) +
  theme_minimal() +
  NULL
permute_plot


plot_name <- file.path("results/variant_table_plots/Null_distribution_from_permutations.png")
ggsave(plot_name, permute_plot, width = 6, height = 4)
