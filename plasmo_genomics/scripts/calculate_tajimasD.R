
# Script to calculate Tajima's D across genome using VCF file
# Output includes: Tajima's D per 5-SNP window, and a plot
# Usage: Rscript calculate_tajimasD.R /path/to/input.vcf /output/dir

# Load libraries
suppressPackageStartupMessages({
  library(vcfR)
  library(pegas)
})

# -------------------- CONFIG --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(" Usage: Rscript calculate_tajimasD.R input.vcf output_directory")
}
vcf_path <- args[1]
out_dir <- args[2]

snps_per_window <- 5  # adjust for low-SNP genomes

# -------------------- FUNCTIONS --------------------
filter_vcf <- function(vcf_obj, max_missing = 0.5) {
  # Filter for biallelic SNPs and low missingness
  vcf_biallelic <- vcf_obj[is.biallelic(vcf_obj), ]
  gt <- extract.gt(vcf_biallelic, element = "GT", as.numeric = FALSE)
  missingness <- apply(gt, 1, function(x) mean(is.na(x) | x == "./."))
  vcf_biallelic[missingness <= max_missing, ]
}

calculate_tajima_sliding <- function(gen_matrix, positions, snps_per_window) {
  n_windows <- ceiling(ncol(gen_matrix) / snps_per_window)
  tajima_vals <- numeric(n_windows)
  midpoints <- numeric(n_windows)
  
  for (i in seq_len(n_windows)) {
    idx <- ((i - 1) * snps_per_window + 1):min(i * snps_per_window, ncol(gen_matrix))
    
    if (length(idx) >= 3) {
      gen_win <- gen_matrix[, idx, drop = FALSE]
      tajima_vals[i] <- tryCatch(tajima.test(gen_win)$D, error = function(e) NA)
      midpoints[i] <- mean(positions[idx])
    } else {
      tajima_vals[i] <- NA
      midpoints[i] <- NA
    }
  }
  
  na.omit(data.frame(midpoint = midpoints, tajima_d = tajima_vals))
}

# -------------------- ANALYSIS --------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_path)

cat("Filtering VCF...\n")
vcf_filtered <- filter_vcf(vcf)

cat("Converting to DNAbin...\n")
gen <- vcfR2DNAbin(vcf_filtered)
positions <- as.numeric(getFIX(vcf_filtered)[, "POS"])

# Ensure genotype matrix and positions align
stopifnot(ncol(gen) == length(positions))

cat("Calculating Tajima's D in", snps_per_window, "SNP windows...\n")
tajima_df <- calculate_tajima_sliding(gen, positions, snps_per_window)

# Save results
out_file <- file.path(out_dir, paste0("tajimas_d_", snps_per_window, "snp_windows.tsv"))
write.table(tajima_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Tajima's D values saved to:", out_file, "\n")

# -------------------- PLOT --------------------
png(file.path(out_dir, "tajimas_d_plot.png"), width = 800, height = 600)
plot(tajima_df$midpoint, tajima_df$tajima_d,
     pch = 20, col = "steelblue",
     xlab = "Genomic Position", ylab = "Tajima's D",
     main = paste("Tajima's D (", snps_per_window, "-SNP Windows)", sep = ""))
abline(h = 0, col = "gray50", lty = 2)
dev.off()

pdf(file.path(out_dir, "tajimas_d_plot.pdf"), width = 8, height = 6)
plot(tajima_df$midpoint, tajima_df$tajima_d,
     pch = 20, col = "steelblue",
     xlab = "Genomic Position", ylab = "Tajima's D",
     main = paste("Tajima's D (", snps_per_window, "-SNP Windows)", sep = ""))
abline(h = 0, col = "gray50", lty = 2)
dev.off()

cat("Tajima's D plots saved.\n")
