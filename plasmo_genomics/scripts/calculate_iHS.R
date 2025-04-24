

# Script to compute and plot iHS using the rehh package
# Usage: Rscript calculate_iHS.R /path/to/rehh_output

# Load libraries
suppressPackageStartupMessages({
  library(rehh)
  library(ggplot2)
})

# ----------------- CONFIG --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript calculate_iHS.R /path/to/rehh_output")
}

base_dir <- args[1]
chrom_dirs <- list.dirs(base_dir, recursive = FALSE)

# ---------------------------------------------
# Function to transpose hap file to rehh format
transpose_hap_file <- function(hap_path, out_path) {
  hap_data <- read.table(hap_path, header = FALSE, stringsAsFactors = FALSE)
  transposed <- t(hap_data)
  hap_ids <- paste0("hap", seq_len(nrow(transposed)))
  transposed_df <- cbind(ID = hap_ids, as.data.frame(transposed))
  write.table(transposed_df, file = out_path, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  message("Transposed hap saved to: ", out_path)
}

# ----------------- iHS CALC -------------------
haplohh_list <- list()

for (chrom_path in chrom_dirs) {
  chrom_name <- basename(chrom_path)
  hap_file <- file.path(chrom_path, paste0(chrom_name, ".hap"))
  map_file <- file.path(chrom_path, paste0(chrom_name, ".map"))
  hap_fixed <- file.path(chrom_path, paste0(chrom_name, ".transposed.hap"))
  
  if (file.exists(hap_file) && file.exists(map_file)) {
    transpose_hap_file(hap_file, hap_fixed)
    cat("Processing", chrom_name, "\n")
    
    haplohh_obj <- data2haplohh(hap_file = hap_fixed,
                                map_file = map_file,
                                allele_coding = "01",
                                min_maf = 0.05,
                                chr.name = chrom_name,
                                verbose = FALSE)
    
    scan_result <- scan_hh(haplohh_obj)
    haplohh_list[[chrom_name]] <- scan_result
  } else {
    warning(" Missing hap/map for ", chrom_name)
  }
}

# ----------------- Combine iHS -------------------
all_ihs_list <- lapply(names(haplohh_list), function(chr) {
  ihs <- ihh2ihs(haplohh_list[[chr]])$ihs
  ihs_df <- as.data.frame(ihs)
  ihs_df$CHR <- chr
  return(ihs_df)
})

combined_ihs <- do.call(rbind, all_ihs_list)
write.table(combined_ihs, file = file.path(base_dir, "combined_iHS_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Combined iHS written to file.\n")

# ----------------- Manhattan Plot -------------------
ihs_data <- na.omit(read.table(file.path(base_dir, "combined_iHS_results.txt"),
                               header = TRUE, stringsAsFactors = FALSE))
names(ihs_data)[names(ihs_data) == "POSITION"] <- "BP"
names(ihs_data)[names(ihs_data) == "IHS"] <- "iHS_raw"
ihs_data$abs_iHS <- abs(ihs_data$iHS_raw)
ihs_data$CHR <- as.factor(ihs_data$CHR)

# Cumulative position
ihs_data <- ihs_data[order(ihs_data$CHR, ihs_data$BP), ]
ihs_data$BPcum <- NA
offset <- 0
tick_pos <- c()

for (chr in unique(ihs_data$CHR)) {
  chr_idx <- ihs_data$CHR == chr
  ihs_data$BPcum[chr_idx] <- ihs_data$BP[chr_idx] + offset
  tick_pos <- c(tick_pos, mean(ihs_data$BPcum[chr_idx]))
  offset <- max(ihs_data$BPcum[chr_idx])
}

# Identify top 1% outliers
cutoff <- quantile(ihs_data$abs_iHS, 0.99)
top_hits <- subset(ihs_data, abs_iHS >= cutoff)

# Plot
manhattan_plot <- ggplot(ihs_data, aes(x = BPcum, y = abs_iHS, color = CHR)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text(data = top_hits, aes(label = BP), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), length.out = length(unique(ihs_data$CHR)))) +
  scale_x_continuous(label = unique(ihs_data$CHR), breaks = tick_pos) +
  labs(x = "Chromosome", y = "|iHS|", title = "Manhattan Plot of |iHS| Scores (Top 1% Labeled)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Save
ggsave(file.path(base_dir, "iHS_manhattan_plot.pdf"), manhattan_plot, width = 8, height = 6)
ggsave(file.path(base_dir, "iHS_manhattan_plot.png"), manhattan_plot, width = 8, height = 6)
cat("Manhattan plot saved to PDF and PNG.\n")
