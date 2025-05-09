# Usage: Rscript fst_plot.R <fst_clean_file> <chr_list> <outlier_output_bed> <annotated_genes_input> <prefix>

args <- commandArgs(trailingOnly = TRUE)
fst_file <- args[1]
chr_list_file <- args[2]
bed_outfile <- args[3]
gene_file <- args[4]
prefix <- args[5]

library(dplyr)
library(ggplot2)
library(ggrepel)
library(fuzzyjoin)
library(stringr)

# Load FST data
fst_data <- read.delim(fst_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(fst_data) <- c("chrom", "midPos", "Nsites", "weighted_fst")

# Filter chromosomes and clean data
chr_list <- readLines(chr_list_file)
fst_data <- fst_data %>%
  filter(chrom %in% chr_list, !is.na(weighted_fst), is.finite(weighted_fst)) %>%
  mutate(chrom = factor(chrom, levels = unique(chrom)))

# Cumulative positions
chrom_info <- fst_data %>%
  group_by(chrom) %>%
  summarize(chr_len = max(midPos), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len)

fst_data <- fst_data %>%
  left_join(chrom_info, by = "chrom") %>%
  mutate(pos = midPos + tot)

# Axis labels
axis_labels <- fst_data %>%
  group_by(chrom) %>%
  summarize(center = (min(pos) + max(pos)) / 2)

# Identify top 1% outliers
fst_threshold <- quantile(fst_data$weighted_fst, 0.99, na.rm = TRUE)
fst_outliers <- fst_data %>%
  filter(weighted_fst >= fst_threshold) %>%
  mutate(start = midPos - 5000, end = midPos + 5000)

# Save BED
write.table(
  fst_outliers[, c("chrom", "start", "end")],
  bed_outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Annotate genes
annotated_genes <- read.delim(gene_file, sep = "\t", stringsAsFactors = FALSE)
annotated_genes <- annotated_genes %>%
  mutate(gene_id = str_replace(gene_id, "ID=gene-", "")) %>%
  select(chrom = chr, gff_start, gff_end, gene_id) %>%
  rename(start = gff_start, end = gff_end)

fst_outliers_annot <- fst_outliers %>%
  group_by(chrom) %>%
  interval_inner_join(
    annotated_genes %>% group_by(chrom),
    by = c("start", "end")
  ) %>%
  ungroup() %>%
  transmute(chrom = chrom.x, midPos, weighted_fst, pos, gene_id) %>%
  distinct()

# Plot Manhattan
ggplot(fst_data, aes(x = pos, y = weighted_fst, color = chrom)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = fst_threshold, color = "red", linetype = "dashed") +
  geom_point(data = fst_outliers_annot, color = "red", size = 2) +
  geom_text_repel(
    data = fst_outliers_annot %>% filter(!is.na(gene_id)),
    aes(label = gene_id), size = 3, max.overlaps = 10
  ) +
  scale_x_continuous(breaks = axis_labels$center, labels = axis_labels$chrom) +
  scale_color_manual(values = rep(c("steelblue", "orange"), length(unique(fst_data$chrom)))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, size = 8)
  ) +
  labs(
    title = paste0("FST Manhattan Plot with Annotated Top 1% Genes (", prefix, ")"),
    x = "Chromosome", y = "Weighted FST"
  )

