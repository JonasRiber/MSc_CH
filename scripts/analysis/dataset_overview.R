
library(tidyverse)
library(data.table)


df <- read.delim("C:\\Users\\au683396\\Documents\\sample_manifest.txt")



#remove doners
df <- df |> filter(patientID_biobank != "na")
df$sample_timepoint_days <- as.numeric(df$sample_timepoint_days)



#finding last sample
df <- df |> 
  group_by(patientID_cluster) |> 
  mutate(last_sample = max(sample_timepoint_days)) |> 
  ungroup()

#reorder based on last sample
df$patientID_cluster <- factor(df$patientID_cluster, levels = df |> 
                          distinct(patientID_cluster, last_sample) |> 
                          arrange(desc(last_sample)) |> 
                          pull(patientID_cluster))


#convert days to years
df <- df |> mutate(years = as.numeric(sample_timepoint_days)/365)

# Plot
ggplot(df, aes(x = sample_timepoint_days, y = patientID_cluster)) +
  geom_point(size = 2, 
             color = "black",
             fill = "skyblue",
             shape = 21,
             stroke = 0.6,
             na.rm = TRUE) +
  geom_line(aes(group = patientID_cluster), color = "gray70", linetype = "dashed") +
  scale_y_discrete(name = "Patients",
                   limits = rev(levels(df$patientID_cluster)))+
  labs(
    title = "Timeline of Sample Collections per Patient",
    x = "Days") +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(-5, 1550)) +
  NULL


# mean sample count
sample_counts <- df |> group_by(patientID_cluster) |> summarise(n = n()) 
mean(sample_counts$n)

# median sample count
median(sample_counts$n)

# total sample count
sum(sample_counts$n)




df_detect <- read.csv("C:\\Users\\au683396\\Documents\\2023_feb_Detection_result.csv", sep = ";")



#remove doners
df_detect <- df_detect |> filter(biobankID != "na")
df_detect$sample_timepoint_days_sinse_OP <- as.numeric(df_detect$sample_timepoint_days_sinse_OP)

df_detect$biobankID <- as.character(df_detect$biobankID)


#finding last sample
df_detect <- df_detect |> 
  group_by(biobankID) |> 
  mutate(last_sample = max(sample_timepoint_days_sinse_OP)) |> 
  ungroup()

#reorder based on last sample
df_detect$biobankID <- factor(df_detect$biobankID, levels = df_detect |> 
                                 distinct(biobankID, last_sample) |> 
                                 arrange(desc(last_sample)) |> 
                                 pull(biobankID))

#convert days to years
df_detect <- df_detect |> mutate(years = as.numeric(sample_timepoint_days_sinse_OP)/365)


# Plot
ggplot(df_detect, aes(x = sample_timepoint_days_sinse_OP, y = biobankID)) +
  geom_point(aes(fill = sample_cat), size = 2, 
             color = "black",
             #fill = "skyblue",
             shape = 21,
             stroke = 0.6) +
  geom_line(aes(group = biobankID), color = "gray70", linetype = "dashed") +
  scale_y_discrete(name = "Patients",
                   limits = rev(levels(df_detect$biobankID)))+
  labs(
    title = "Timeline of Sample Collections per Patient",
    x = "Days",
    fill = "Sample \ncategory") +
  #scale_color_manual(
  #  values = c("act" = "#1F77B4", "post-act" = "#FF7F0E", "post-op" = "#2CA02C", "pre-op" = "#D62728", "surveillance" = "#9467BD"),
  #  labels = c("act" = "Act", "post-act" = "Post-act", "post-op" = "Post-op", "pre-op" = "Pre-op", "surveillance" = "Surveillance")
  #)+
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(-5, 1100)) +
  NULL






#sample counts
nsamp_detect <- df_detect |> group_by(biobankID) |> summarise(n = n())
mean(nsamp_detect$n)  #mean
median(nsamp_detect$n)#median
sum(nsamp_detect$n)   #total
sum(df_detect$TUMOR_FRACTION_spec==0) #non ctDNA detected


