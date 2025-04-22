library(ggplot2)
library(reshape2)
library(scales)

# Read sample names
sample_names <- read.table("pruned_variants_filtered.fam", header = FALSE, stringsAsFactors = FALSE)[,1]

# Read ADMIXTURE Q matrix (adjust the K value as needed)
admix_data <- read.table("pruned_variants_filtered.5.Q", header = FALSE)

# Ensure length match
stopifnot(nrow(admix_data) == length(sample_names))

# Merge sample names and extract MD codes
admix_data <- cbind(Sample = sample_names, admix_data)
admix_data$Sample <- gsub(".*-(MD[0-9]+).*", "\\1", admix_data$Sample)

# Rename clusters
colnames(admix_data) <- c("Sample", paste0("Cluster_", 1:(ncol(admix_data) - 1)))

# Reshape
admix_long <- melt(admix_data, id.vars = "Sample", variable.name = "Cluster", value.name = "Ancestry")

# Plot
ggplot(admix_long, aes(x = Sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  labs(
    x = "Isolates",
    y = "Ancestry (%)",
    title = "Population Structure Analysis of Pm Isolates (ADMIXTURE K=5)"
  )
