
library(tidyverse)
library(data.table)


#test patient: PATEINT_ID
patient_id <- "PATIENT_ID"
sample_id <- "SAMPLE_ID"


gformat <- function(x, gscale="Mb"){
  # genomic scale: bp, kb, Mb, Gb
  base_gscale <- c("bp"=0, "kb"=3, "mb"=6, "gb"=9)
  if( !(tolower(gscale) %in% names(base_gscale)) ){
    stop(paste("gscale not valid. Must be one of:", paste(names(base_gscale), collapse = ", ")))}
  y <- x/10^base_gscale[tolower(gscale)]
  z <- paste(format(y, digits=3, big.mark=",", scientific=FALSE), gscale)
  # z <- paste(sprintf("%.3g", y), gscale) # "%.1f"
  return(z)
}




### single sample SAMPLE_ID ###
df1_single <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", sample_id,"_cfdna_depth.tsv"), sep = "\t", header = FALSE)
setDT(df1_single)

df1_single_mms <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", sample_id,"_cfdna_mismatches.tsv"), sep = "\t", header = FALSE)
setDT(df1_single_mms)


df1_mean <- round(mean(df1_single$V5), 2)
df1_n <- nrow(df1_single)
df1_tot_mm <- sum(df1_single_mms$V8) 
df1_vaf <- NA

plot1 <- df1_single |> ggplot(aes(x = V5)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=1)+
  coord_cartesian(ylim = c(0, 0.1))+
  labs(y = "") +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot1


text_box_data <- data.frame(
  x = c(1, 2, 3,4),
  y = 0.5,
  label = c(paste0(gformat(df1_n), "*"), 
            paste0(gformat(df1_tot_mm), "*"),
            paste0(df1_mean, "X"),
            paste0(df1_vaf)),
  title = c("Total sites", "Total mms", "Mean depth", "Mean VAF")
)

#* = since the depth is collected af having merge all the mismatch files together.
# the same number of sites are found in the single and merge samples, in the single sample
# some of these sites are not mismatched, hence the total number of sites is bigger than the actual 
# mismatches

annotation_plot1 <- ggplot(text_box_data, aes(x = x, y = y)) +
  geom_tile(width = 1, height = 1, fill = "#B0F26D", color = "black") +
  geom_text(aes(label = label), size = 3, lineheight = 1) +
  geom_text(aes(y=y+0.65, label = title), size = 3) +
  theme_void() +  # clean background
  theme(plot.margin = margin(10, 10, 10, 10))
annotation_plot1




### Merged depth ###
df2_merged_depth <- read.table(file =  paste0("C:\\Users\\au683396\\Documents\\", patient_id,"_merged_mnv_calls_depth.tsv"), sep = "\t", header = FALSE)

df2_mean <- round(mean(df2_merged_depth$V9), 2)
df2_n <- nrow(df2_merged_depth)
df2_tot_mm <- sum(df2_merged_depth$V9)
df2_vaf <- round(mean(df2_merged_depth$V10),3)

plot2 <- df2_merged_depth |> ggplot(aes(x = V9)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50)+
  labs(y = "Proportion of sites") +
  coord_cartesian(xlim = c(100, 600)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),)
plot2


text_box_data <- data.frame(
  x = c(1, 2, 3, 4),
  y = 0.5,
  label = c(paste0(gformat(df2_n)), 
            paste0(gformat(df2_tot_mm)),
            paste0(df2_mean, "X"),
            paste0(df2_vaf))
)

annotation_plot2 <- ggplot(text_box_data, aes(x = x, y = y)) +
  geom_tile(width = 1, height = 1, fill = "#B0F26D", color = "black") +
  geom_text(aes(label = label), size = 3, lineheight = 1) +
  theme_void() +  # clean background
  theme(plot.margin = margin(10, 10, 10, 10))
annotation_plot2






### Filtering of the data ### 
df4_somatic_calls <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", patient_id, "_somatic_calls_w_depth.tsv"), sep = "\t", header = FALSE)
names(df4_somatic_calls) <- c("chr", "start", "end", "ref", "alt", "Mtype", ".", "n_alt", "depth","vaf", "som_threshold")
setDT(df4_somatic_calls)

upper_depth <- mean(df4_somatic_calls$depth) + 3*sd(df4_somatic_calls$depth)

df4_somatic_calls[n_alt>3 & depth>100 & depth<upper_depth] |> 
  ggplot(aes(x = depth)) +
  geom_histogram(aes(y=after_stat(count / sum(count))), bins=50)



### basic filtering ### 
df4_somatic_calls <- df4_somatic_calls |> filter(depth>100 & depth<upper_depth & n_alt>3)

df3_mean <- round(mean(df4_somatic_calls$depth),3)
df3_n <- nrow(df4_somatic_calls)
df3_tot_mm <- sum(df4_somatic_calls$n_alt)
df3_vaf <- round(mean(df4_somatic_calls$vaf),3)

plot3 <- df4_somatic_calls[n_alt>3 & depth>100 & depth<upper_depth] |> 
  ggplot(aes(x = depth)) +
  geom_histogram(aes(y=after_stat(count / sum(count))), bins=50) +
  labs(x = "",
       y = "")+
  #caption = "100<depth<mean+3sd\nn_alt>3\nsomatic removed") +
  coord_cartesian(xlim = c(100, 600)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot3


text_box_data <- data.frame(
  x = c(1, 2, 3, 4),
  y = 0.5,
  label = c(paste0(gformat(df3_n, gscale = "kb")), 
            paste0(gformat(df3_tot_mm)),
            paste0(df3_mean, "X"),
            paste0(df3_vaf))
)

annotation_plot3 <- ggplot(text_box_data, aes(x = x, y = y)) +
  geom_tile(width = 1, height = 1, fill = "#B0F26D", color = "black", ) +
  geom_text(aes(label = label), size = 3, lineheight = 1) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))
annotation_plot3




#somatic threshold
df4_somatic_calls <- df4_somatic_calls |> mutate(group = ifelse(test = n_alt>som_threshold, yes = "above", no = "below"))
df4_somatic_calls <- df4_somatic_calls |> filter(group == "below")

df4_mean <- round(mean(df4_somatic_calls$depth),3)
df4_n <- nrow(df4_somatic_calls)
df4_tot_mm <- sum(df4_somatic_calls$n_alt)
df4_vaf <- round(mean(df4_somatic_calls$vaf),3)


plot4 <- df4_somatic_calls[group == "below"] |> 
  ggplot(aes(x = depth)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50)+
  labs(x = "Depth",
       y = "")+
       #caption = "100<depth<mean+3sd\nn_alt>3\nsomatic removed") +
  coord_cartesian(xlim = c(100, 600)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 0),
        axis.title.y = element_blank())
plot4


text_box_data <- data.frame(
  x = c(1, 2, 3, 4),
  y = 0.5,
  label = c(paste0(df4_n), 
            paste0(gformat(df4_tot_mm, gscale = "kb")),
            paste0(df4_mean, "X"),
            paste0(df4_vaf))
  )

annotation_plot4 <- ggplot(text_box_data, aes(x = x, y = y)) +
  geom_tile(width = 1, height = 1, fill = "#B0F26D", color = "black", ) +
  geom_text(aes(label = label), size = 3, lineheight = 1) +
  theme_void() +  # clean background
  theme(plot.margin = margin(10, 10, 10, 10))
annotation_plot4







### full plot ###
library(patchwork)
row_title <- function(label) {
  ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = label, size = 3,
             angle = 90) +
    theme_void()
}


row1 <- row_title("Single sample") + (plot1 + annotation_plot1) + 
  plot_layout(widths = c(0.1, 1))  # ← was 0.1, now 0.05

row2 <- row_title("Aggregated \nsamples") + (plot2 + annotation_plot2) + 
  plot_layout(widths = c(0.1, 1))

row3 <- row_title("Static filter") + (plot3 + annotation_plot3) + 
  plot_layout(widths = c(0.1, 1))

row4 <- row_title("Somatic filter") + (plot4 + annotation_plot4) + 
  plot_layout(widths = c(0.1, 1))

ytitle <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Proportion of sites", size = 4,
           angle = 90) +
  theme_void() +
  plot_layout(widths = c(0.01,0.5))

full_plot <- ytitle + (row1 / row2 / row3 / row4) 
full_plot








df4_somatic_calls[group == "below" & n_alt>3] |> 
  ggplot(aes(x = depth, fill = chr)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50)+
  labs(title = "Chromosomal presense in the bimodal peaks", 
       x = "Depth",
       y = "Proportion of sites",
       caption = "100<depth<mean+3sd\nn_alt>3\nsomatic removed",
       fill = "Chromosome") +
  annotate(geom = "text", x = 500, y = 0.095, 
           label = paste0("Total sites: ", df4_n, 
                          "\nMean depth: ", df4_mean
           )) +
  coord_cartesian(xlim = c(100, 550), ylim = c(0, 0.1)) +
  theme_minimal() 


### investigating bimodalality
#split the distributions
df <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", patient_id, "_somatic_calls_w_depth_w_buffy.tsv"), sep = "\t", header = FALSE)
names(df) <- c("chr", "start", "end", "ref", "alt", "Mtype", ".", "n_alt", "depth","vaf", "som_threshold")
setDT(df)
df <- df |> mutate(group = ifelse(test = n_alt>som_threshold, yes = "above", no = "below"))
df <- df |> filter(chr != "chrX")


middle <- 260
df <- df |> mutate(modality = ifelse(test = depth>middle, yes = "above", no = "below"))
df <- df |> filter(depth>100 & depth<(mean(df$depth)+3*sd(df$depth)))


# only removing germline
df[group == "below"] |> 
  ggplot(aes(x = depth)) +
  geom_histogram(binwidth=1)+
  labs(title = "3.filtered calls", 
       x = "depth",
       y = "proportion of sites") +
  theme_minimal() +
  NULL

# Showcasing chromomsome
df[group == "below" & n_alt>3] |> 
  ggplot(aes(x = depth, fill = chr)) +
  geom_vline(xintercept = middle, linetype = "dashed")+
  geom_histogram(bins=50)+
  labs(title = "Chromosomes across depth distribution", 
       x = "Depth",
       y = "Proportion of sites") +
  theme_minimal() +
  NULL


# lower peak presence of chr
df[modality=="below"] |> group_by(chr) |> summarise(n_positions = n()) |> 
  ggplot(aes(x = chr, y = n_positions, fill = chr)) +
  geom_col() +
  labs(title = "Chromosomal presence in lower peak",
       x = "Chromosome",
       y = "n") +
  theme_minimal()+
  theme(legend.position = "none") + 
  NULL

df[modality=="above"] |> group_by(chr) |> summarise(n_positions = n()) 


#
ggplot() +
  geom_histogram(data=df, aes(x = depth, y = after_stat(count/sum(count))), fill = "red", alpha = 0.5) +
  geom_histogram(data=df[group=="below" & n_alt>3], aes(x = depth, y = after_stat(count/sum(count))),
                 fill = "blue", alhpa = 0.2) +
  NULL

df[group=="below" & n_alt>7] |>
  ggplot(aes(x = vaf))+
  geom_histogram(aes(y = after_stat(count / sum(count)))) + 
  facet_wrap(~modality) +
  NULL


df[group=="below" & n_alt>3] |> ggplot(aes(x = depth)) +
  geom_histogram()




#chrX placement
dfX <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", patient_id, "_somatic_calls_w_depth_no_dups_w_X.tsv"), 
                  sep = "\t", header = FALSE)


names(dfX) <- c("chr", "start", "end", "ref", "alt", "Mtype", ".", "n_alt", "depth","vaf", "som_threshold")
setDT(dfX)

dfX <- dfX |> mutate(group = ifelse(test = n_alt>som_threshold, yes = "above", no = "below"))
dfX <- dfX |> mutate(chr_state = ifelse(test = chr=="chrX", yes = "chrX", no = "Not chrX"))


dfX <- dfX |> filter(depth<(mean(depth)+3*sd(depth)) & depth>100)


dfX_mean <- round(mean(dfX$depth), 2)
dfX_n <- nrow(dfX)

chrX_mean<- mean(dfX[chr=="chrX"]$depth)

pt1 <- dfX[group=="below"] |> 
  ggplot(aes(x = depth, fill = chr_state)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50)+
  geom_vline(xintercept = chrX_mean, linetype = "dashed", color = "red") +
  labs(title = "Only depth filtering - no nr. mm fitler", 
       x = "Depth",
       y = "Proportion of sites",
       caption = "100<depth<mean+3sd\nsomatic removed") +
  annotate(geom = "text", x = 450, y = 0.09, 
           label = paste0("Total sites: ", dfX_n, 
                          "\nMean depth: ", dfX_mean
          )) +
  coord_cartesian(xlim = c(100, 600), ylim = c(0, 0.1)) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  NULL
pt1


dfX_mean <- round(mean(dfX[group=="below" & n_alt>3]$depth), 2)
dfX_n <- nrow(dfX[group=="below" & n_alt>3])

pt2 <- dfX[group=="below" & n_alt>3] |> 
  ggplot(aes(x = depth, fill = chr_state)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50)+
  geom_vline(xintercept = chrX_mean, linetype = "dashed", color = "red") +
  labs(title = "With nr. mm>3 filter", 
       x = "Depth",
       y = "Proportion of sites",
       caption = "100<depth<mean+3sd\nn_alt>3\nsomatic removed") +
  annotate(geom = "text", x = 450, y = 0.09, 
           label = paste0("Total sites: ", dfX_n, 
                          "\nMean depth: ", dfX_mean
           )) +
  coord_cartesian(xlim = c(100, 600), ylim = c(0, 0.1)) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  NULL
pt2

library(patchwork)
pt1 + pt2



dfX <- dfX |> mutate(cutoff = ifelse(test = n_alt>3, yes = "above", no = "below"))

dfX[group=="below" & chr!="chrX" & cutoff == "below"] |> 
  ggplot(aes(x = depth, fill = cutoff)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins=50) +
  labs(title = "sites under mm cutoff (n_alt>3)",
       y = "proportion of sites")+
  NULL

dfX[group=="below" & chr!="chrX" & cutoff == "above"] |> 
  ggplot(aes(x = depth)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 bins=50, fill = "cornflowerblue") +
  labs(title = "sites over mm cutoff (n_alt>3)",
       y = "proportion of sites")+
  NULL





#annotated
df_annotated <- read.table(file = paste0("C:\\Users\\au683396\\Documents\\", patient_id, "_annotated.tsv"), sep = "\t", header = T)
setDT(df_annotated)

middle <- 260

df_annotated <- df_annotated |> mutate(modality = ifelse(test = Total_count>middle, yes = "above", no = "below"))

df_annotated |> group_by(modality, Func) |> 
  summarise(n=n(), .groups = "drop") |>
  group_by(modality) |> 
  mutate(prop = n / sum(n)) |> 
  ggplot(aes(x = Func, y = prop, group = modality, fill = modality)) +
  geom_col(position = "dodge") +
  labs(title = "Ratio of region type within \nthe lower and upper peaks",
       x = "Region type", 
       y = "Proportion") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("above" = "#1f77b4", "below" = "#ff7f0e"),
                    labels = c("above" = "Higher Peak", "below" = "Lower Peak")) +
  NULL



genes_in_both <- df_annotated |> 
  distinct(Gene, modality) |> 
  count(Gene) |> 
  filter(n<2)


genes_not_in_both <- genes_in_both$Gene


variants_group_counts <- df_annotated[AA_change!="."] |> 
  distinct(AA_change, modality) |> 
  count(AA_change) |> 
  filter(n<2)


variant_not_in_both <- variants_group_counts$AA_change


#all annotated calls
df_annotated |> ggplot(aes(x = Total_count, fill=Gene)) +
  geom_histogram() +
  theme(legend.position = "none")

#only exonic
df_annotated[AA_change!="."] |> 
  ggplot(aes(x = Total_count, fill = Gene)) +
  geom_histogram()


                                    




