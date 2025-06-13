

#### Author: Jonas Riber Jørgensen
#### date: 11/04/2025
#### description: R-script to process and visualize the filtered calls

### initialization
library(tidyverse)
library(data.table)

setwd("/path/to/wokring/directory")

if(!interactive()) pdf(NULL) # suppress Rplots.pdf file


### setup args
args <- commandArgs(trailingOnly = TRUE)

#check if a file path was provided
if (length(args) < 2) {
  stop("Usage: Rscript visualize_calls.R <input_file.tsv> <output_directory>")
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


setDT(data)


file_name <- tools::file_path_sans_ext(basename(file_path))
id <- print(str_extract(string = file_name, pattern = "ID_PATTERN"))


#get meta data
metadata <- read.delim("data/Metadata.tsv", header=TRUE, sep="\t")

#metadata <- read.delim("C:\\Users\\au683396\\Documents\\Metadata.tsv", header=TRUE, sep="\t")
sample_count <- metadata$sample_count[metadata$patientID == id]









### plot 6 - Using the selected filter values and adding drivers
library(openxlsx)
driver_list <- read.xlsx(xlsxFile = "data/blood_bld-2022-018825-mmc1.xlsx")

#local
#driver_list <- read.xlsx(xlsxFile = "C:\\Users\\au683396\\Documents\\blood_bld-2022-018825-mmc1.xlsx")

#main drivers from https://www.nature.com/articles/s41467-022-31878-0
main_drivers <- data.frame(Gene = c("DNMT3A", "TET2", "ASXL1", "TP54", "JAK2", "SF3B1"))


#mark somatic or germ-snp
data <- data |> mutate(group = ifelse(Alt_count < Somatic_threshold, "somatic", "germline"))

pvalue <- 1e-10


#check if any gene is in the driver list
data <- data |> 
  mutate(
    category = sapply(strsplit(Gene, ","), function(x) {
      if (any(x %in% main_drivers$Gene)) {
        return("main_driver")
      } else if (any(x %in% driver_list$Gene)) {
        return("driver")
      } else {
        return("passenger")
      }
    })
  )

mean_depth <- mean(data$Total_count) 





### somatic vs germline on the set filter values
somatic_count <- data |> group_by(group) |> summarise(n=n())
somatic_count$n[1]

plot <- data |>
  ggplot(aes(x = VAF, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 70) +
  geom_text(aes(x = 0.4, y = 800 * 1.1, label = paste0("n=", somatic_count$n[1])), 
            vjust = 2, size = 3, color = "firebrick") + 
  geom_text(aes(x = 0.1, y = 700 * 1.1, label = paste0("n=", somatic_count$n[2])), 
            vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Zoomed to X:0-0.45 and Y:0-900\n", 
                        "Mean depth: ", mean_depth, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>100 & n_alt>3")) +
  #coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 900)) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_6.1somatic_vs_germ_VAF_distribution_filtered.png"))
ggsave(plot_name, plot)


### isolating somatic with driver genes - showing "the bump"
plot <- data[group=="somatic"] |>
  ggplot(aes(x = VAF, fill = category)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
  geom_text(aes(x = 0.1, y = 100000 * 1.1, label = paste0("n=", somatic_count$n[2])), 
            vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", mean_depth, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>100 & n_alt>3")) +
  coord_cartesian(xlim = c(0, 0.45), ylim = c(0,50)) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_6.2somatic_VAF_distribution_w_drivers_filtered.png"))
ggsave(plot_name, plot)



### depth distribution
plot <- data[group=="somatic"] |> ggplot(aes(x = Total_count))+
  geom_histogram() +
  #geom_vline(xintercept = 260) +
  labs(title = "Depth distribution",
       x = "Depth",
       y = "n",
       caption = paste0("Mean depth: ", mean_depth, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>100 & n_alt>3")) +
  NULL
plot


plot_name <- file.path(output_dir, paste0(id, "_7depth_distribution_filtered.png"))
ggsave(plot_name, plot)


summary(data$Total_count)






### bimodal depth investigation
# old data before exclussion

#file_path <- "C:\\Users\\au683396\\Documents\\_annotated_old.tsv"
#output_dir <- getwd()
#data_old <- read.delim(file_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#above
#setDT(data)
#data[Total_count>260] |> group_by(Func) |> summarise(n = n())
#nrow(data[Total_count>260])

#below
#data[Total_count<260] |> group_by(Func) |> summarise(n = n()) |> 
#  ggplot(aes(x = Func, y = n)) +
#  geom_col()
#nrow(data[Total_count<260])

#depth_threshold <- 260
#data <- data |> mutate(state=ifelse(Total_count>depth_threshold, "above", "below"))
#data_old <- data_old |> mutate(state=ifelse(Total_count>depth_threshold, "above", "below"))
#data_old <- data_old |> mutate(group = ifelse(Alt_count < Somatic_threshold, "somatic", "germline"))


#data_summary <- data |> 
#  group_by(state, Func) |> 
#  summarise(n = n(), .groups = "drop") |> 
#  group_by(state) |> 
#  mutate(prop = n / sum(n))


#data_summary |> 
#  ggplot(aes(x = Func, y = prop, fill = state)) +
#  geom_col(position = "dodge")
#
#data_summary <- data |> 
#  group_by(state, group) |> 
#  summarise(n = n(), .groups = "drop") |> 
#  group_by(state) |> 
#  mutate(prop = n / sum(n))

#data_summary |> 
#  ggplot(aes(x = group, y = prop, fill = state)) +
#  geom_col(position = "dodge")



#setDT(data_old)

#library(patchwork)
#no_germ <- data_old |> ggplot(aes(x = Total_count)) +
#  geom_histogram() +
#  geom_vline(xintercept = depth_threshold) +
#  coord_cartesian(ylim = c(0, 1000)) +
#  labs(title = "No germline filter",
#       x = "Depth")

#som_test <- data_old[group == "somatic"] |> ggplot(aes(x = Total_count)) +
#  geom_histogram() +
#  geom_vline(xintercept = depth_threshold) +
#  coord_cartesian(ylim = c(0, 1000)) +
#  labs(title = "snp removel - binom test",
#       x = "Depth")

#snp_exclussion <- data |> ggplot(aes(x = Total_count)) +
#  geom_histogram() +
#  geom_vline(xintercept = depth_threshold) +
#  coord_cartesian(ylim = c(0, 1000)) +
#  labs(title = "snp sxclussion",
#       x = "Depth")


#both <- data[group == "somatic"] |>  
#  ggplot(aes(x = Total_count)) +
#  geom_histogram() +
#  geom_vline(xintercept = depth_threshold) +
#  coord_cartesian(ylim = c(0, 1000)) +
#  labs(title = "both",
#       x = "Depth")

# no_germ + som_test + snp_exclussion + both + labs(caption = "Panel 1: No filtering done for germline
#                                                   Panel 2: removel via the binomial test  and somatic threshold
#                                                   Panel 3: Removal via the called snps
#                                                   Panel 4: combination of panel 2 and 3")

 
 



#data[group == "somatic" & Alt_count>7] |>  
#   ggplot(aes(x = Total_count)) +
#   geom_histogram() +
#   geom_vline(xintercept = depth_threshold) +
#   coord_cartesian(ylim = c(0, 150)) +
#   labs(title = "",
#        x = "")
 
 
#data[group=="somatic"] |> 
#  ggplot(aes(x = Start, y = Total_count)) +
#  geom_point()+
#  #coord_cartesian(ylim = c(0, 150)) +
#  NULL




#above
#setDT(data_old)
#data_old[Total_count>260] |> group_by(Func) |> summarise(n = n())
#nrow(data_old[Total_count>260])

#below
#data_old[Total_count<260] |> group_by(Func) |> summarise(n = n())
#nrow(data_old[Total_count<260])
#data_summary_old <- data_old |> 
#  group_by(state, Func) |> 
#  summarise(n = n(), .groups = "drop") |> 
#  group_by(state) |> 
#  mutate(prop = n / sum(n))

#data_summary_old |> 
#  ggplot(aes(x = Func, y = prop, fill = state)) +
#  geom_col(position = "dodge")


#data_summary_old <- data_old |> 
#  group_by(state, group) |> 
#  summarise(n = n(), .groups = "drop") |> 
#  group_by(state) |> 
#  mutate(prop = n / sum(n))

#data_summary_old |> 
#  ggplot(aes(x = group, y = prop, fill = state)) +
#  geom_col(position = "dodge")


#df <- left_join(data_summary, data_summary_old, by = c("state", "Func"))
#df_long <- df %>%
#  pivot_longer(cols = c(prop.x, prop.y),
#               names_to = "group",
#               values_to = "prop") %>%
#  mutate(group = recode(group, 
#                        "prop.x" = "After",
#                        "prop.y" = "Before"))


#df_long |> ggplot(aes(x = Func, y = prop, group = group, fill = group)) +
#  geom_col(position = "dodge")

#x <- data |> mutate(depth = Total_count-Cons_count)





### Gene representation
#gene_data <- data |> 
#  filter(group == "somatic") |> 
#  group_by(category, Gene) |> 
#  separate_rows(Gene, sep = ",") |>
#  count(Gene, name = "n") 

#gene_data |> 
#  arrange(desc(n)) |> 
#  slice_head(n = 15) |> 
#  ggplot(aes(x = fct_reorder(Gene, n, .desc = FALSE), y = n, fill=category)) +
#  geom_col() +
#  labs(x = "Gene", y = "Count") +
#  theme_minimal() +
#  coord_flip() +
#  NULL




#driver_data <- data[category!="passenger"] |> 
#  filter(group == "somatic") |> 
#  group_by(category, Gene) |> 
#  separate_rows(Gene, sep = ",") |>
#  count(Gene, name = "n") 

#driver_data |>
#  arrange(desc(n)) |> 
#  slice_head(n = 15) |> 
#  ggplot(aes(x = fct_reorder(Gene, n, .desc = FALSE), y = n, fill=category)) +
#  geom_col() +
#  labs(x = "Gene", y = "Count") +
#  theme_minimal() +
#  coord_flip() +
#  NULL






