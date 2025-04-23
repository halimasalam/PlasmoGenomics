# ADMIXTURE Plotting Script for Population Structure Analysis

# Load required libraries
library(ggplot2)     # For plotting
library(reshape2)    # For data reshaping (wide to long format)
library(scales)      # For formatting y-axis labels as percentages

# Path to the .fam file (used to get sample names)
fam_file <- "analysis_results/admixture/files/pruned_variants_filtered.fam"

# Path to the Q-matrix file from ADMIXTURE (change K as needed)
# Example: If using K=5, the Q file will be `*.5.Q`
q_file <- "analysis_results/admixture/files/pruned_variants_filtered.5.Q"

# Path to the Output plot filename
output_plot_file <- "analysis_results/admixture/plot/admix_barplot_K5.png"    

# Read only the first column of the .fam file (sample identifiers)
sample_names <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)[, 1]

# Load ADMIXTURE Q-matrix data
admix_data <- read.table(q_file, header = FALSE)

# Ensure that the number of rows in Q-matrix matches the number of samples
stopifnot(nrow(admix_data) == length(sample_names))

# Merge sample names with Q-matrix and process labels
# Add sample names as a new column
admix_data <- cbind(Sample = sample_names, admix_data)

# Extract MD identifiers (e.g., "MD9", "MD12") from sample names
admix_data$Sample <- gsub(".*-(MD[0-9]+).*", "\\1", admix_data$Sample)

# Rename ancestry cluster columns as Cluster_1, Cluster_2, ..., Cluster_K
num_clusters <- ncol(admix_data) - 1
colnames(admix_data) <- c("Sample", paste0("Cluster_", 1:num_clusters))

# Reshape data for plotting
# Convert from wide format (one row per sample) to long format
admix_long <- melt(admix_data, id.vars = "Sample", 
                   variable.name = "Cluster", value.name = "Ancestry")

# Generate ADMIXTURE bar plot
admix_plot <- ggplot(admix_long, aes(x = Sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +  # Stacked bar plot
  scale_y_continuous(labels = percent_format()) +  # Show ancestry as percentages
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Rotate sample labels
    legend.position = "none"  # Hide legend (optional)
  ) +
  labs(
    x = "Isolates",
    y = "Ancestry (%)",
    title = paste("Population Structure Analysis (ADMIXTURE K =", num_clusters, ")")
  )


# Display the plot
print(admix_plot)

# Save as PNG
ggsave(output_plot_file, plot = admix_plot, width = 10, height = 6, dpi = 300)

