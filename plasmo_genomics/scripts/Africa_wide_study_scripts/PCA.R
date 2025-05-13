# ------------------- Load Required Libraries -------------------
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# ------------------- FUNCTION DEFINITIONS -------------------

# Load genotype matrix and format it
load_genotype_matrix <- function(gt_file, sample_file) {
  gt <- read.table(gt_file, header = FALSE, stringsAsFactors = FALSE)
  sample_ids <- readLines(sample_file)
  colnames(gt) <- c("SNP", sample_ids)
  
  gt_t <- as.data.frame(t(gt[, -1]))
  colnames(gt_t) <- gt$SNP
  rownames(gt_t) <- sample_ids
  
  return(gt_t)
}

# Convert genotype characters to numeric
convert_to_numeric <- function(gt_df) {
  convert_genotype <- function(genotype) {
    if (genotype == "./.") return(NA)
    return(as.numeric(genotype))
  }
  numeric_df <- as.data.frame(lapply(gt_df, function(x) sapply(x, convert_genotype)))
  rownames(numeric_df) <- rownames(gt_df)
  return(numeric_df)
}

# Filter SNPs and impute missing values
filter_and_impute <- function(gt_numeric, remove_outliers = TRUE, outlier_threshold = 200) {
  var_filter <- apply(gt_numeric, 2, var, na.rm = TRUE)
  filtered <- gt_numeric[, !is.na(var_filter) & var_filter != 0]
  
  imputed <- as.data.frame(apply(filtered, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  if (remove_outliers) {
    pca_temp <- prcomp(imputed, center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca_temp$x)
    outliers <- which(abs(pca_df$PC1) > outlier_threshold)
    imputed <- imputed[-outliers, ]
    
    # Refilter after removing outliers
    var_check <- apply(imputed, 2, var, na.rm = TRUE)
    imputed <- imputed[, !is.na(var_check) & var_check != 0]
    
    # Re-impute just in case
    imputed <- as.data.frame(apply(imputed, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    }))
  }
  
  return(imputed)
}

# Perform PCA
run_pca <- function(gt_matrix) {
  pca_result <- prcomp(gt_matrix, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SampleID <- rownames(pca_df)
  percent_var <- round(100 * summary(pca_result)$importance[2, 1:2], 2)
  return(list(pca = pca_result, coords = pca_df, var = percent_var))
}

# Merge metadata and PCA results
merge_with_metadata <- function(pca_coords, metadata_file) {
  metadata <- read_tsv(metadata_file, show_col_types = FALSE)
  merged <- left_join(pca_coords, metadata, by = "SampleID")
  return(merged)
}

# Plot PCA by Country
plot_pca_by_country <- function(pca_df, percent_var) {
  country_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(pca_df$Country)))
  
  ggplot(pca_df, aes(x = PC1, y = PC2, fill = Country)) +
    geom_point(shape = 21, size = 3, alpha = 0.8) +
    scale_fill_manual(values = country_colors) +
    theme_minimal() +
    coord_cartesian(xlim = c(-120, 120), ylim = c(-350, 350)) +
    labs(
      title = "PCA of Plasmodium vivax Samples by Country",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      fill = "Country"
    )
}

# Plot PCA by Region
plot_pca_by_region <- function(pca_df, title_prefix = "Plasmodium malariae") {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Region)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste("PCA of", title_prefix, "Isolates by Region"),
      x = paste0("PC1 (", round(100 * var(pca_df$PC1) / sum(apply(pca_df[, c("PC1", "PC2")], 2, var)), 1), "%)"),
      y = paste0("PC2 (", round(100 * var(pca_df$PC2) / sum(apply(pca_df[, c("PC1", "PC2")], 2, var)), 1), "%)")
    ) +
    scale_color_brewer(palette = "Set2") +
    coord_cartesian(xlim = c(-60, 60), ylim = c(-60, 60))
}

# ------------------- EXECUTION -------------------

# File paths (update these as needed)
gt_file <- "path/to/gt_matrix.tsv"
sample_file <- "path/to/sample_names.txt"
metadata_file <- "path/to/metadata.tsv"  ##(sample/country/region)

# Pipeline execution
gt_raw <- load_genotype_matrix(gt_file, sample_file)
gt_numeric <- convert_to_numeric(gt_raw)
gt_clean <- filter_and_impute(gt_numeric, remove_outliers = TRUE)

pca_results <- run_pca(gt_clean)
pca_coords <- pca_results$coords
percent_var <- pca_results$var

# Merge with metadata
pca_merged <- merge_with_metadata(pca_coords, metadata_file)

# Plot by country
plot_pca_by_country(pca_merged, percent_var)

# Plot by region (alternative)
plot_pca_by_region(pca_merged, title_prefix = "Plasmodium vivax")
