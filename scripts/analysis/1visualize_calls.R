

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



setnames(data, old = c("X.CHROM", "START","END", "REF", "ALT", "MTYPE", "N_CON", "N_ALT", "N_ALL", 
                       "VAF", "X"), 
         new = c("chr", "start", "end", "ref", "alt", "Mtype", "n_con", "n_alt", "depth", 
                 "VAF", "somatic_threshold"))
setDT(data)



file_name <- tools::file_path_sans_ext(basename(file_path))
id <- print(str_extract(string = file_name, pattern = "ID_PATTERN"))


#get meta data
metadata <- read.delim("data/Metadata.tsv", header=TRUE, sep="\t")

#metadata <- read.delim("C:\\Users\\au683396\\Documents\\Metadata.tsv", header=TRUE, sep="\t")
sample_count <- metadata$sample_count[metadata$patientID == id]


### data summary
head(data)

summary(data$depth)

merged_mean_depth <- round(mean(data$depth), 2)

nrow(data)
nrow(data[depth>5])
nrow(data[depth>20])
nrow(data[depth>50])
nrow(data[depth>100])
nrow(data[n_alt>1])

nrow(data[depth>20&n_alt>1])



#############
### Plots ###
#############

### plot1.1: distribution of VAF
plotx <- data |> ggplot(aes(x = VAF)) +
  geom_histogram(fill = "steelblue") +
  labs(title = "VAF distribution",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth , 
                        "\nSample count: ", sample_count, 
                        "\nNo filtering")) +
  theme_minimal() +
  NULL
plotx

plot_name <- file.path(output_dir, paste0(id, "_1.1VAF_distribution.png"))
ggsave(plot_name, plot)


### plot1.2: distribution of VAF
merged_mean_depth_20 <- mean(data[depth>20]$depth)
plot <- data[depth>20] |> ggplot(aes(x = VAF)) +
  geom_histogram(fill = "steelblue") +
  labs(title = "VAF distribution - filtered",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth_20 , 
                        "\nSample count: ", sample_count, 
                        "\nFiltering: depth>20")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_1.2VAF_distribution_filtered.png"))
ggsave(plot_name, plotx)



### plot1.3: distribution of VAF
max_depth <- round(mean(data$depth)+3*sd(data$depth), 2)

merged_mean_depth_100 <- round(mean(data[depth>100 & n_alt>3]$depth), 2)
#data[VAF < VAF_upp & depth > 100 & n_alt>3 & depth<max_depth]


ploty <- data[depth>100&n_alt>3] |> ggplot(aes(x = VAF)) +
  geom_histogram(fill = "steelblue") +
  labs(title = "VAF distribution - filtered",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth_100 , 
                        "\nSample count: ", sample_count, 
                        "\nFiltering: ", "depth>100",
                        "\nmin nr alternative: 3")) +
  theme_minimal()+ 
  #xlim()
  NULL
ploty

library(patchwork)

plotx / ploty
plot_name <- file.path(output_dir, paste0(id, "_1.3VAF_distribution_filtered2.png"))
ggsave(plot_name, ploty)


### plot1.4: distribution of VAF - zoomed to 0.02-0.4
VAF_upp <- 0.4

plot <- data[VAF < VAF_upp & depth > 100 & n_alt>3 & depth<max_depth] |> 
  ggplot(aes(x = VAF)) +
  geom_histogram(fill = "steelblue") +
  labs(title = "VAF distribution - zoomed",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth_100 , 
                        "\nSample count: ", sample_count, 
                        "\nFiltering: VAF<", VAF_upp, "and depth>100",
                        "\nFiltering: ", max_depth,">depth>100",
                        "\nmin nr alternative: 3")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_1.4VAF_distribution_filtered2_zoomed.png"))
ggsave(plot_name, plot)






### plot2: density distribution of VAF 
plot <- data[depth>100] |> ggplot(aes(x = VAF)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(title = "Density Plot of VAF",
       x = "VAF (alt/depth)",
       y = "Density",
       caption = paste0("Mean depth: ", merged_mean_depth_100 , 
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>100")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_2VAF_density_distribution.png"))
ggsave(plot_name, plot)




### plot3.1: scatterplot
pvalue <- 1e-9
data <- data |> mutate(group = ifelse(n_alt < somatic_threshold, "Somatic", "Germline"))

data[depth>100 & n_alt>3 & group == "somatic"]

plot <- data |> ggplot(aes(x = VAF, y = depth, group = group, colour = group)) +
  geom_point(alpha = 0.6) +
  labs(title = "VAF vs. Read Depth",
       x = "Variant Allele Frequency",
       y = "Total Read Depth",
       color = "H0 count",
       caption = paste0(
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nNo filtering")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_3.1VAF_vs_depth.png"))
ggsave(plot_name, plot)



### plot3.2: distribution of lower than H0
plot <- data[group == "somatic"] |> 
  ggplot(aes(x = VAF)) +
  geom_histogram(fill = "steelblue") +
  labs(title = "VAF distribution",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Only variants with lower n_alt count than somatic_threshold\n",
                        "Mean depth: ", merged_mean_depth, 
                        "\nqbinom pvalue: ", pvalue, 
                        "\nSample count: ", sample_count, 
                        "\nNo filtering")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_3.2distribution_lower_H0.png"))
ggsave(plot_name, plot)


### plot 3.3: distribution of higher than H0
plot <- data[group == "germline"] |> 
  ggplot(aes(x = VAF)) +
  geom_histogram(fill ="firebrick") +
  labs(title = "VAF distribution",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Only variants with higher somatic count than H0_count\n",
                        "Mean depth: ", merged_mean_depth, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nNo filtering")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_3.3distribution_higher_H0.png"))
ggsave(plot_name, plot)




data[depth > 100 & n_alt>3 & group == "Somatic" & depth<max_depth] |> 
  ggplot(aes(x = VAF)) +
  geom_histogram() +
  NULL




### plot 4.1: VAF distribution - grouping colored
somatic_count <- data[depth > 20] |> group_by(group) |> summarise(n=n())
somatic_count$n[1]

plot <- data[depth > 20] |>
  ggplot(aes(x = VAF, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  geom_text(aes(x = 0.8, y = 5000 * 1.1, label = paste0("n=", somatic_count$n[1])), 
            vjust = 2, size = 3, color = "firebrick") + 
  geom_text(aes(x = 0.2, y = 3000 * 1.1, label = paste0("n=", somatic_count$n[2])), 
            vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth_20, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>20")) +
  NULL
plot

plot_name <- file.path(output_dir, paste0(id, "_4.1VAF_distribution_qbinom_grouping.png"))
ggsave(plot_name, plot)


### plot 4.2: VAF distribution - grouped colored and zoomed
plot <- data[depth > 20] |>
  ggplot(aes(x = VAF, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  geom_text(aes(x = 0.4, y = 800 * 1.1, label = paste0("n=", somatic_count$n[1])), 
            vjust = 2, size = 3, color = "firebrick") + 
  geom_text(aes(x = 0.1, y = 700 * 1.1, label = paste0("n=", somatic_count$n[2])), 
            vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Zoomed to X:0-0.5 and Y:0-900\n", 
                        "Mean depth: ", merged_mean_depth_20, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>20")) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 900)) +
  NULL
plot


plot_name <- file.path(output_dir, paste0(id, "_4.2VAF_distribution_qbinom_grouping_zoomed.png"))
ggsave(plot_name, plot)


### plot 4.3: VAF distribution - grouping colored
somatic_count <- data[depth > 100 & depth < max_depth & n_alt > 3] |> group_by(group) |> summarise(n=n())
somatic_count$n[1]

plotx <- data[depth > 100 & depth < max_depth & n_alt > 3] |>
  ggplot(aes(x = VAF, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
  #geom_text(aes(x = 0.8, y = 5000 * 1.1, label = paste0("n=", somatic_count$n[1])), 
  #          vjust = 2, size = 3, color = "firebrick") + 
  #geom_text(aes(x = 0.2, y = 3000 * 1.1, label = paste0("n=", somatic_count$n[2])), 
  #          vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Mean depth: ", merged_mean_depth_100, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: ", max_depth,">depth>100",
                        "\nmin nr of alternative: 3")) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  NULL
plotx

plot_name <- file.path(output_dir, paste0(id, "_4.3VAF_distribution_qbinom_grouping_filter2.png"))
ggsave(plot_name, plotx)


### plot 4.4: VAF distribution - grouped colored and zoomed
ploty <- data[depth > 100  & depth < max_depth & n_alt > 3] |>
  ggplot(aes(x = VAF, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
  geom_text(aes(x = 0.4, y = 800 * 1.1, label = paste0("n=", somatic_count$n[1])), 
            vjust = 2, size = 3, color = "firebrick") + 
  geom_text(aes(x = 0.1, y = 700 * 1.1, label = paste0("n=", somatic_count$n[2])), 
            vjust = 2, size = 3, color = "cornflowerblue") + 
  labs(title = "",
       x = "VAF (alt/depth)",
       y = "n",
       caption = paste0("Zoomed to X:0-0.5 and Y:0-900\n", 
                        "Mean depth: ", merged_mean_depth_100, 
                        "\nqbinom pvalue: ", pvalue,
                        "\nSample count: ", sample_count, 
                        "\nFilter: depth>100")) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 100)) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  NULL
ploty


plot_name <- file.path(output_dir, paste0(id, "_4.4VAF_distribution_qbinom_grouping_zoomed_filter2.png"))
ggsave(plot_name, ploty)


plotx / ploty


### plot 5 - depth overview for different filter levels
library(patchwork)


median_depth <- median(data$depth)
mean_depth <- mean(data$depth)
upper_cutoff <- mean(data$depth) + 3 * sd(data$depth)


depth_plot <- data |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff*1.1, y = nrow(data)/20,
           label = round(upper_cutoff, 2),,
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff, linetype = "dashed", color = "red")+
  labs(subtitle = "No filter",
       caption = paste0("Mean: ", mean_depth, 
                        "\nMedian:", median_depth,
                        "\nSample count: ", sample_count))



depth_data1 <- data[depth>20]
median_depth1 <- median(depth_data1$depth)
mean_depth1 <- mean(depth_data1$depth)
upper_cutoff1 <- mean(depth_data1$depth) + 3 * sd(depth_data1$depth)


depth_plot1 <- depth_data1 |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff1*1.1, y = nrow(data[depth>20])/20,
           label = round(upper_cutoff1, 2),,
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff1, linetype = "dashed", color = "red")+
  labs(subtitle = "Filter: depth>20",
       caption = paste0("Mean: ", mean_depth1, 
                        "\nMedian:", median_depth1,
                        "\nSample count: ", sample_count))



depth_data2 <- data[n_alt>1]
median_depth2 <- median(depth_data2$depth)
mean_depth2 <- mean(depth_data2$depth)
upper_cutoff2 <- mean(depth_data2$depth) + 3 * sd(depth_data2$depth)

depth_plot2 <- depth_data2 |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff2*1.1, y = nrow(data)/50,
           label = round(upper_cutoff2, 2),,
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff2, linetype = "dashed", color = "red")+
  labs(subtitle = "Filter: n_alt>1",
       caption = paste0("Mean: ", mean_depth2, 
                        "\nMedian:", median_depth2,
                        "\nSample count: ", sample_count))


depth_data3 <- data[depth>20 & n_alt>1]
median_depth3 <- median(depth_data3$depth)
mean_depth3 <- mean(depth_data3$depth)
upper_cutoff3 <- mean(depth_data3$depth) + 3 * sd(depth_data3$depth)


depth_plot3 <- depth_data3 |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff3*1.1, y = nrow(data[depth>20])/50,
           label = round(upper_cutoff3, 2),,
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff3, linetype = "dashed", color = "red")+
  labs(subtitle = "Filter: depth>20 & n_alt>1",
       caption = paste0("Mean: ", mean_depth3, 
                        "\nMedian:", median_depth3,
                        "\nSample count: ", sample_count))


depth_data4 <- data[depth>50 & n_alt>1]
median_depth4 <- median(depth_data4$depth)
mean_depth4 <- mean(depth_data4$depth)
upper_cutoff4 <- mean(depth_data4$depth) + 3 * sd(depth_data4$depth)

depth_plot4 <- depth_data4 |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff4*1.1, y = nrow(data[depth>50])/50,
           label = round(upper_cutoff4, 2),,
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff4, linetype = "dashed", color = "red")+
  labs(subtitle = "Filter: depth>50 & n_alt>1",
       caption = paste0("Mean: ", mean_depth4, 
                        "\nMedian:", median_depth4,
                        "\nSample count: ", sample_count))


depth_data5 <- data[depth>100 & n_alt>3]
median_depth5 <- median(depth_data5$depth)
mean_depth5 <- round(mean(depth_data5$depth), 2)
upper_cutoff5 <- mean(depth_data5$depth) + 3 * sd(depth_data5$depth)

depth_plot5 <- depth_data5 |> ggplot(aes(x = depth))+
  geom_histogram() +
  annotate("text", x = upper_cutoff5*1.1, y = nrow(data[depth>100])/50,
           label = round(upper_cutoff5, 2),
           size=3, color="red") +
  geom_vline(xintercept = upper_cutoff5, linetype = "dashed", color = "red")+
  labs(subtitle = "Filter: depth>100 & n_alt>3",
       caption = paste0("Mean: ", mean_depth5, 
                        "\nMedian:", median_depth5,
                        "\nSample count: ", sample_count))
depth_plot5

# combine plots
plot <- depth_plot + depth_plot1 + depth_plot2 + depth_plot3 + depth_plot4 + depth_plot5
plot

plot_name <- file.path(output_dir, paste0(id, "_5depth_overview_w_different_filters.png"))
ggsave(plot_name, plot)

plot_name <- file.path(output_dir, paste0(id, "_5.1depth_100rd_3nalt.png"))
ggsave(plot_name, depth_plot5)




### plot 6 - Using the selected filter values and adding drivers
#library(openxlsx)
#driver_list <- read.xlsx(xlsxFile = "data/blood_bld-2022-018825-mmc1.xlsx")

#local
#driver_list <- read.xlsx(xlsxFile = "C:\\Users\\au683396\\Documents\\blood_bld-2022-018825-mmc1.xlsx")

#main drivers - https://www.nature.com/articles/s41467-022-31878-0
#main_drivers <- data.frame(Gene = c("DNMT3A", "TET2", "ASXL1", "TP54", "JAK2", "SF3B1"))



#filter data - no max_depth
#filtered_data <- data[depth>100]
#print(paste0("Full = ", nrow(data),"    Filtered = ", nrow(filtered_data)))

#check if any gene is in the driver list
#filtered_data <- filtered_data |> 
#  mutate(
#    category = sapply(strsplit(gene, ","), function(x) {
#      if (any(x %in% main_drivers$Gene)) {
#        return("main_driver")
#      } else if (any(x %in% driver_list$Gene)) {
#        return("driver")
#      } else {
#        return("passenger")
#      }
#    })
#  )

#filtered_mean_depth <- mean(filtered_data$depth) 




### somatic vs germline on the set filter values
#somatic_count <- filtered_data |> group_by(group) |> summarise(n=n())
#somatic_count$n[1]

#plot <- filtered_data |>
#  ggplot(aes(x = VAF, fill = group)) +
#  geom_histogram(position = "identity", alpha = 0.7, bins = 70) +
#  geom_text(aes(x = 0.4, y = 800 * 1.1, label = paste0("n=", somatic_count$n[1])), 
#            vjust = 2, size = 3, color = "firebrick") + 
#  geom_text(aes(x = 0.1, y = 700 * 1.1, label = paste0("n=", somatic_count$n[2])), 
#            vjust = 2, size = 3, color = "cornflowerblue") + 
#  labs(title = "",
#       x = "VAF (alt/depth)",
#       y = "n",
#       caption = paste0("Zoomed to X:0-0.45 and Y:0-900\n", 
#                        "Mean depth: ", filtered_mean_depth, 
#                        "\nqbinom pvalue: ", pvalue,
#                        "\nSample count: ", sample_count, 
#                        "\nFilter: depth>50")) +
#  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 900)) +
#  NULL
#plot

#plot_name <- file.path(output_dir, paste0(id, "_6.1somatic_vs_germ_VAF_distribution_filtered.png"))
#ggsave(plot_name, plot)


### isolating somatic with driver genes - showing "the bump"
#plot <- filtered_data[group=="somatic"] |>
#  ggplot(aes(x = VAF, fill = category)) +
#  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
#  geom_text(aes(x = 0.1, y = 100000 * 1.1, label = paste0("n=", somatic_count$n[2])), 
#            vjust = 2, size = 3, color = "cornflowerblue") + 
#  labs(title = "",
#       x = "VAF (alt/depth)",
#       y = "n",
#       caption = paste0("Mean depth: ", filtered_mean_depth, 
#                        "\nqbinom pvalue: ", pvalue,
#                        "\nSample count: ", sample_count, 
#                        "\nFilter: depth > 50")) +
#  coord_cartesian(xlim = c(0, 0.45), ylim = c(0, 800)) +
#  NULL
#plot

#plot_name <- file.path(output_dir, paste0(id, "_6.2somatic_VAF_distribution_w_drivers.png"))
#ggsave(plot_name, plot)


### filter max_depth aswell
#filtered_data <- filtered_data[depth<mean(depth)+3*sd(depth)]
#print(paste0("Full = ", nrow(data),"    Filtered = ", nrow(filtered_data)))
#filtered_mean_depth <- mean(filtered_data$depth) 

#somatic_count <- filtered_data |> group_by(group) |> summarise(n=n())
#somatic_count$n[1]

#plot <- filtered_data[group=="somatic"] |>
#  ggplot(aes(x = VAF, fill = category)) +
#  geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
#  geom_text(aes(x = 0.1, y = 100000 * 1.1, label = paste0("n=", somatic_count$n[2])), 
#            vjust = 2, size = 3, color = "cornflowerblue") + 
#  labs(title = "",
#       x = "VAF (alt/depth)",
#       y = "n",
#       caption = paste0("Mean depth: ", filtered_mean_depth, 
#                        "\nqbinom pvalue: ", pvalue,
#                        "\nSample count: ", sample_count, 
#                        "\nFilter:mean(depth)+3*sd(depth) > depth > 50")) +
#  coord_cartesian(xlim = c(0, 0.45), ylim = c(0, 800)) +
#  NULL
#plot

#plot_name <- file.path(output_dir, paste0(id, "_6.3somatic_VAF_distribution_w_drivers_and_max_depth.png"))
#ggsave(plot_name, plot)

#note: the bump now dissapears, investigating this in IGV shows that some positions that
# the blacklist aimed to remove are still present (areas with excessive coverage)
# additionally multi-allelic positions seem to be placed in the section aswell







####################
### experimental ###
####################



#head(data[depth > 100 & VAF>0.8 & VAF<0.98 & chr!="chrX"], 10)


#data[end>1805000 & end<1805500]

#data[VAF>1]







#gene_data <- data |> 
#  mutate(
#    category = sapply(strsplit(gene, ","), function(x) {
#      if (any(x %in% main_drivers$Gene)) {
#        return("main_driver")
#      } else if (any(x %in% driver_list$Gene)) {
#        return("driver")
#      } else {
#        return("passenger")
#      }
#    })
#  )

#gene_mean_depth <- mean(gene_data$depth) 


#summary(gene_data$depth)



#driver_data <- gene_data[category=="driver" | category=="main_driver"] |> 
#  filter(depth<500 & depth>20 & n_alt>1 & VAF>0.02 & group == "somatic") |> 
#  group_by(gene) |> 
#  separate_rows(gene, sep = ",") |>
#  count(gene, name = "n") 
  
#driver_data |> 
#  ggplot(aes(x = fct_reorder(gene, n, .desc = FALSE), y = n)) +
#  geom_col() +
#  labs(x = "Gene", y = "Count") +
#  theme_minimal() +
#  coord_flip() 




#x <- filtered_data[depth<mean(depth)+3*sd(depth) & depth>200 & group=="somatic" & VAF>0.1 & VAF<0.4]

#driver_data <- filtered_data[category=="driver" | category=="main_driver"] |> 
#  filter(VAF>0.1 & VAF<0.4 & group=="somatic")

### driver genes
#filtered_data[category=="driver" | category=="main_driver"] |>
#  filter(group=="somatic")

#plot <- filtered_data[category=="driver" | category=="main_driver"] |>
#  ggplot(aes(x = VAF, fill = category)) +
#  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  #geom_text(aes(x = 0.4, y = 800 * 1.1, label = paste0("n=", somatic_count$n[1])), 
  #          vjust = 2, size = 3, color = "firebrick") + 
  #geom_text(aes(x = 0.1, y = 700 * 1.1, label = paste0("n=", somatic_count$n[2])), 
  #          vjust = 2, size = 3, color = "cornflowerblue") + 
#  labs(title = "",
#       x = "VAF (alt/depth)",
#       y = "n",
#       caption = paste0("Zoomed to X:0-0.5 and Y:0-900\n", 
#                        "Mean depth: ", merged_mean_depth_100, 
#                        "\nqbinom pvalue: ", pvalue,
#                        "\nSample count: ", sample_count, 
#                        "\nFilter: depth>100 & n_alt>1")) +
#  #coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 900)) +
#  NULL
#plot



### VAF vs depth - adding the drivers
#driver_snps <- filtered_data[category=="driver"]

#plot <- filtered_data |> 
#  ggplot() +
#  geom_point(aes(x = VAF, y = depth, group=group, colour = group), alpha = 0.6) +
#  geom_point(mapping = aes(x = VAF, y = depth), data = driver_snps, colour = "yellow") +
#  labs(title = "VAF vs. Read Depth",
#       x = "Variant Allele Frequency",
#       y = "Total Read Depth",
#       color = "H0 count",
#       caption = paste0("Mean depth: ", merged_mean_depth_100 ,
#                        "\nqbinom pvalue: ", pvalue,
#                        "\nSample count: ", sample_count, 
#                        "\nFilter: depth>50 & n_alt>1")) +
#  NULL
#plot



#nrow(data[depth>50 & n_alt>1 & group=="somatic"])


#old <- 3990009
#new <- 217402
  
#(old-new)/old*100


#as.data.table(data)[depth>=20 & n_alt>1, .(vaf=mean(VAF), n=.N), by=.(group, exonic_func)]

#exo_data <- data[depth>20 & n_alt>1 & func=="exonic" & group=="somatic"]
#mean(exo_data$VAF)




#0.1–1 10−2



