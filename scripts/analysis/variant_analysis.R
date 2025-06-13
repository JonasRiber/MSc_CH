


library(tidyverse)
library(data.table)

df <- read.table(file = "C:\\Users\\au683396\\Documents\\variant_table_mms_depth.tsv", sep = "\t", header = TRUE)

# get patient cols
pt <- grep("PT_ID", names(df), value = TRUE)

# function to convert "X | Y" -> X / Y
convert_to_vaf <- function(x) {
  sapply(strsplit(x, "\\|"), function(parts) {
    num <- as.numeric(trimws(parts[1]))
    den <- as.numeric(trimws(parts[2]))
    if (!is.na(num) && !is.na(den) && den != 0) {
      return(num / den)
    } else {
      return(NA)
    }
  })
}

#apply the conversion only to selected columns
df[pt] <- lapply(df[pt], convert_to_vaf)

setDT(df)
#pivot longer
df_vaf <- df |> melt(measure.vars=pt, value.name = "vaf", variable.name="pt")



pt_count <- length(unique(df_vaf$pt))
df_vaf <- df_vaf  |>  mutate(gene = str_extract(AA_change, "^[^:]+"))
df_vaf <- df_vaf |> mutate(odd_SETBP = ifelse(AA_change == "SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T", yes = "odd_var", no = "not_odd"))

#--- plots ---#
df_vaf |> filter(vaf > 0) |> 
  ggplot(aes(x=vaf, fill = odd_SETBP))+ 
  geom_histogram(aes(y = after_stat(count / sum(count)))) +
  facet_wrap(~pt) +
  theme_bw() +
  NULL


df_vaf |> filter(vaf>0) |> 
  ggplot(aes(x = vaf)) +
  geom_histogram(aes( y = after_stat(count / sum(count)))) +
  facet_wrap(~gene)+
  theme_bw()+
  labs(title = "VAF Distributions within each Gene",
       x = "VAF",
       y = "Proportion",
       caption = "Entire cohort, filter=vaf>0") +
  NULL



### singular variants
#SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T
#looking at the table itself, it seems like SETBP1 has way more mismatches than other variants
#the two most mismatched variants are from this gene.

SETBP1_counts <- df_vaf |> filter(AA_change == "SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T") |>
  summarise(total = n(), zero_vaf = sum(vaf==0), above_zero = total-zero_vaf)
SETBP1_counts

plot <- df_vaf |> filter(AA_change == "SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T")  |> 
  ggplot(aes(x=vaf)) +
  geom_histogram(fill = "firebrick") +
  labs(title = "SETBP1 vaf distribution",
       x = "VAF",
       y = "n",
       caption = paste0("VAF distribution of a single variant across all patients
       SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T\n",
                        "nr of zero: ", SETBP1_counts[[2]][1], 
                        "\nnr of non-zero: ", SETBP1_counts[[3]][1])) +
  theme_minimal() +
  NULL
plot

plot_name <- file.path("results/variant_table_plots/SETBP1_single_vaf.png")
ggsave(plot_name, plot)



#"SETBP1:NM_001130110.2:exon4:c.A682T:p.T228S"
SETBP1_counts2 <- df_vaf |> filter(AA_change == "SETBP1:NM_001130110.2:exon4:c.A682T:p.T228S") |>
  summarise(total = n(), zero_vaf = sum(vaf==0), above_zero = total-zero_vaf)
SETBP1_counts2


df_vaf |> filter(AA_change == "SETBP1:NM_001130110.2:exon4:c.A682T:p.T228S")  |> 
  ggplot(aes(x=vaf)) +
  geom_histogram(fill = "firebrick") +
  labs(title = "SETBP1 vaf distribution",
       x = "VAF",
       y = "n",
       caption = paste0("VAF distribution of a single variant across all patients
       SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T\n",
                        "nr of zero: ", SETBP1_counts2[[2]][1], 
                        "\nnr of non-zero: ", SETBP1_counts2[[3]][1])) +
  theme_minimal() +
  NULL



#"DNMT3A:NM_001320893.1:exon13:c.G1703T:p.R568L,DNMT3A:NM_001375819.1:exon13:c.G1490T:p.R497L,DNMT3A:NM_153759.3:exon14:c.G1592T:p.R531L,DNMT3A:NM_022552.5:exon18:c.G2159T:p.R720L,DNMT3A:NM_175629.2:exon18:c.G2159T:p.R720L"
SETBP1_counts2 <- df_vaf |> filter(AA_change == "DNMT3A:NM_001320893.1:exon13:c.G1703T:p.R568L,DNMT3A:NM_001375819.1:exon13:c.G1490T:p.R497L,DNMT3A:NM_153759.3:exon14:c.G1592T:p.R531L,DNMT3A:NM_022552.5:exon18:c.G2159T:p.R720L,DNMT3A:NM_175629.2:exon18:c.G2159T:p.R720L")  |> 
  summarise(total = n(), zero_vaf = sum(vaf==0), above_zero = total-zero_vaf)
SETBP1_counts2

df_vaf |> filter(AA_change == "DNMT3A:NM_001320893.1:exon13:c.G1703T:p.R568L,DNMT3A:NM_001375819.1:exon13:c.G1490T:p.R497L,DNMT3A:NM_153759.3:exon14:c.G1592T:p.R531L,DNMT3A:NM_022552.5:exon18:c.G2159T:p.R720L,DNMT3A:NM_175629.2:exon18:c.G2159T:p.R720L")  |> 
  ggplot(aes(x=vaf)) +
  geom_histogram(fill = "firebrick") +
  labs(title = "DNMT3A vaf distribution",
       x = "VAF",
       y = "n",
       caption = paste0("VAF distribution of a single variant across all patients
       DNMT3A:NM_001320893.1:exon13:c.G1703T:p.R568L...\n",
                        "nr of zero: ", SETBP1_counts2[[2]][1], 
                        "\nnr of non-zero: ", SETBP1_counts2[[3]][1])) +
  theme_minimal() +
  NULL





#--- within genes ---#
df_vaf$gene <- sub(":.*", "", df_vaf$AA_change)


genes <- c("SETBP1", "CREBBP", "SETDB1", "DNMT3A", "TET2", "ASXL1")

annotations <- df_vaf |> filter(gene %in% genes & vaf != 0) |> 
  group_by(gene) |> 
  summarise(
    count = n_distinct(pt),
    label1 = paste("Mean VAF:", round(mean(vaf), 2)),
    label2 = paste("Nr patients:", count),
    x = Inf,
    y = Inf
  )


plot <- df_vaf |> filter(gene %in% genes & vaf != 0)|> 
  ggplot(aes(x=vaf, fill = gene)) +
  geom_histogram(aes(y = after_stat(count/sum(count)))) +
  geom_text(data = annotations, aes(x = x, y = y, label = label1), inherit.aes = FALSE,
            size = 3.5, hjust = 1.05, vjust = 1.2, family = "mono") +
  geom_text(data = annotations, aes(x = x, y = y, label = label2), inherit.aes = FALSE,
            size = 3.5, hjust = 1.05, vjust = 3, family = "mono") +
  labs(title = "VAF distribution of individual genes",
       x = "VAF",
       y = "Proportion",
       caption = paste0("VAF distribution of entire gene across all ", pt_count, " patients")) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~gene) +
  coord_cartesian(clip = "off") +
  NULL
plot

plot_name <- file.path("results/variant_table_plots/VAF_of_top_genes.png")
ggsave(plot_name, plot)



### SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T
# removal
df_SETBP <- df_vaf[AA_change == "SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T"] |> filter(vaf!= 0)
summary(df_SETBP$vaf)
length(unique(df_SETBP$pt))








#--- pvalue table ---#
df_pval <- read.table(file = "C:\\Users\\au683396\\Documents\\variant_table_var_ER_pvalues_10req.tsv", header = TRUE)
setDT(df_pval)

df_pval |> filter(AA_change != "SETBP1:NM_001130110.2:exon4:c.G664A:p.A222T")

#load paper ids
id_translation <- read.table("C:\\Users\\au683396\\Documents\\PatientIDs_paper.txt", header = TRUE, sep = "\t", stringsAsFactors = F)


variant_names <- df_pval$AA_change
df_adj <- as.data.frame(apply(df_pval[,..pt,], 2, function(col) p.adjust(col, method = "BH")))

df_adj <- cbind(AA_change = variant_names, df_adj)

sig_counts <- colSums(df_adj[, -1] < 0.1, na.rm = TRUE)

sig_counts <- data.frame(
  patient = names(sig_counts),
  significant_variants = sig_counts
)

setDT(sig_counts)
#sig_counts <- sig_counts |> mutate(pt_id = paste0("pt", row_number()))

sig_counts <- merge(sig_counts, id_translation, by.x = "patient", by.y = "ClusterID", all.x = TRUE)


plot <- sig_counts |> ggplot(aes(x = fct_reorder(PaperID, -significant_variants), y = significant_variants)) +
  geom_col(fill = "firebrick") +
  labs(title = "Significant counts per patient",
       x = "Patients",
       y = "Sig counts",
       caption = "spectrum of how many significant variant each patient
       has before fusing the pvalues. Values have been
       BH adjusted and a threshold of 0.1 was used") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  NULL
plot

plot_name <- file.path("results/variant_table_plots/Significant_variants_per_patient.png")
ggsave(plot_name, plot)



plot <- sig_counts |> group_by(significant_variants) |> 
  summarise(n = n()) |> 
  ggplot(aes(x = significant_variants, y = n)) +
  geom_col(fill = "firebrick")+
  labs(title = "Summation of nr of significant counts",
       x = "Significant variants",
       y = "Number of patients", 
       caption = "Summation of pvalues before fusing them per patient
       adjusted with BH and 0.1 threshold was used") +
  theme_minimal() +
  NULL
plot

plot_name <- file.path("results/variant_table_plots/Significant_variants_per_patient_summation.png")
ggsave(plot_name, plot)


sig_counts |> arrange(-significant_variants)






#--- CH-variants and age ---#
# by number of total variants
variant_count <- df_vaf[vaf != 0] |> group_by(pt) |> summarise(n = n())

# get age
metadata <- read.table(file = "C:\\Users\\au683396\\Documents\\crc__age_sex.txt", sep = "\t", header = TRUE)
metadata$patient_ID <- paste0("PT_ID", metadata$patient_ID)

variant_count <- merge(variant_count, metadata, by.x = "pt", by.y = "patient_ID", all.x = T)
variant_count

#age distribution
summary(variant_count$age)
mean_age <- round(mean(variant_count$age), 2)
variant_count |> ggplot(aes(x = age)) +
  geom_histogram(fill = "cornflowerblue", binwidth = 2) +
  labs(title = "Age distribution of cohort",
       x = "Age",
       y = "n",
       caption = paste0("Mean age ",mean_age)) +
  theme_minimal() +
  NULL

model <- lm(n ~ age, data = variant_count)
slope <- round(coef(model)["age"], 2)
r_squared <- round(summary(model)$r.squared, 2)

variant_count |> ggplot(aes(x = age, y = n)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(subtitle = "Age vs number of different CH variants per patient",
       x = "Age",
       y = "n",
       caption = paste0("Slope: ", slope,
                        "\nR_squared: ", r_squared)) +
  theme_minimal() +
  NULL

variant_count 


# by number of significant variants
sig_counts <- merge(sig_counts, metadata, by.x = "patient", by.y = "patient_ID", all.x = T)
sig_counts

sig_counts |> ggplot(aes(x = age, y = significant_variants)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(title = "Age vs number of different significant variants per patient") +
  NULL









