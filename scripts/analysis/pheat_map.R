

library(tidyverse)
library(data.table)
library(pheatmap)


### based on pvalues
table_path <- "results/variant_table_var_ER_pvalues.tsv"

#local
table_path <- "C:\\Users\\au683396\\Documents\\variant_table_var_ER_pvalues_10req.tsv"

id_translation <- read.table("C:\\Users\\au683396\\Documents\\PatientIDs_paper.txt", header = TRUE, sep = "\t", stringsAsFactors = F)




df <- read.delim(table_path, row.names = 1, check.names = FALSE)
df <- df[, 1:(ncol(df)-3)]


#change patient ids to paperids
oldcols <- colnames(df)
new_cols <- id_translation$PaperID[match(oldcols, id_translation$ClusterID)]
colnames(df) <- new_cols

# add gene
df$gene <- sub(":.*", "", rownames(df))

num_patients <- ncol(df) - 1
patient_cols <- 1:(ncol(df) - 1)  # Assuming 'gene' is last

rownames(df) <- substr(rownames(df), 1, 5)
df$gene <- NULL

library(grid)

p <- pheatmap(as.matrix(df[, !names(df) %in% "gene"]),
         cluster_rows = F,
         cluster_cols = FALSE,
         show_rownames = F,
         show_colnames = T,
         main = "Variant patient table - ptable",
         fontsize = 10,
         silent = TRUE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1, 40), "lines"))))


pushViewport(viewport(layout.pos.col = 1))
grid::grid.text("Variants", x = 0, rot = 90, gp = gpar(fontsize=12))
popViewport()

pushViewport(viewport(layout.pos.col = 2))
grid.draw(p$gtable)
popViewport()

### number of pvals below 1 per variant
below_1_counts <- rowSums(df < 1, na.rm = TRUE)

variant_summary <- data.frame(
  gene = rownames(df),
  num_below_1 = below_1_counts
)

head(variant_summary[order(-variant_summary$num_below_1), ], 10)

sum(variant_summary$num_below_1>1)

### number of pvals below 1 per patient
below_1_counts_cols <- colSums(df < 0.05, na.rm = TRUE)

patient_summary <- data.frame(
  patient = colnames(df),
  num_below_1 = below_1_counts_cols
)

patient_summary[order(-patient_summary$num_below_1), ]



















### based on VAFs
table_path <- "C:\\Users\\au683396\\Documents\\variant_table_mms_depth_10req.tsv"

df_vaf <- read.delim(table_path, row.names = 1, check.names = FALSE)
df_vaf <- df_vaf[, 1:(ncol(df_vaf)-3)]
# get patient cols
pt <- grep("^C0", names(df_vaf), value = TRUE)

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
df_vaf[pt] <- lapply(df_vaf[pt], convert_to_vaf)

setDT(df_vaf)

df_vaf$gene <- sub("_.*", "", rownames(df_vaf))

num_patients <- ncol(df_vaf) - 1
patient_cols <- 1:(ncol(df_vaf) - 1)  # Assuming 'gene' is last
colnames(df_vaf)[patient_cols] <- paste0("pt", seq_along(patient_cols))

rownames(df_vaf) <- substr(rownames(df_vaf), 1, 5)
df_vaf$gene <- NULL


nonzero_breaks <- seq(0.001, 0.16, length.out = 99)
break_list <- c(0, nonzero_breaks)

custom_colors <- c("lightgray",colorRampPalette(c("dodgerblue", "red"))(99))

pheatmap(as.matrix(df_vaf),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = F,
         show_colnames = T,
         color = custom_colors,
         breaks = break_list,
         main = "Variant patient table - VAF")







### How many variants arrive in the table from the variable regions?
#variable regions
#- ZBT33 domain p.9-126, 332-591
#- PPM1D missense exon 5-6
#- CBLB nonsense
#- CBL missense p.381-421
#- ASXL1 and 2 exon 11-12
#- TET2 p1104-1481 and 1843-2002
regions <- data.frame(
  gene = c("ZBT33", "ZBT33", "PPM1D", "CBL", "ASXL1", "ASXL2", "TET2", "TET2"),
  type = c("protein", "protein", "exon", "protein", "exon", "exon", "protein", "protein"),
  start = c(9, 332, 5, 381, 11, 11, 1104, 1843),
  end = c(126, 591, 6, 421, 12, 12, 1481, 2002)
)

get_genes <- function(row) {
  parts <- unlist(strsplit(row, ","))
  genes <- sapply(strsplit(parts, ":"), `[`, 1)
  paste(unique(genes), collapse = ";")
}

get_positions <- function(row, type = "c") {
  parts <- unlist(strsplit(row, ","))
  pattern <- if (type == "c") "c\\.[A-Z0-9]+" else "p\\.[A-Z0-9]+"
  matches <- regmatches(parts, gregexpr(pattern, parts))
  paste(unique(unlist(matches)), collapse = ";")
}


get_exons <- function(row) {
  exon_match <- regmatches(row, regexpr("exon[0-9]+", row))
  exon <- as.numeric(gsub("exon", "", exon_match))
  paste(unique(unlist(exon)), collapse = ";")
}

df$Gene <- sapply(rownames(df), get_genes)
df$cDNA_Change <- sapply(rownames(df), get_positions, type = "c")
df$Protein_change <- sapply(rownames(df), get_positions, type = "p")
df$Exons <- sapply(rownames(df), get_exons)



for (i in 1:nrow(regions)) {
  region <- regions[i, ]
  gene <- region$gene
  type <- region$type
  start <- region$start
  end <- region$end
  
  matches <- rep(FALSE, nrow(df))
  
  for (j in 1:nrow(df)) {
    gene_match <- grepl(gene, df$Gene[j])  # allow semicolon-separated genes
    
    if (!gene_match) next
    
    if (type == "exon") {
      exon_nums <- as.numeric(unlist(strsplit(df$Exons[j], ";")))
      if (any(!is.na(exon_nums) & exon_nums >= start & exon_nums <= end)) {
        matches[j] <- TRUE
      }
    }
    
    if (type == "protein") {
      prot_nums <- as.numeric(gsub("p\\.([A-Z])([0-9]+)[A-Z]", "\\2", unlist(strsplit(df$Protein_change[j], ";"))))
      if (any(!is.na(prot_nums) & prot_nums >= start & prot_nums <= end)) {
        matches[j] <- TRUE
      }
    }
  }
}


