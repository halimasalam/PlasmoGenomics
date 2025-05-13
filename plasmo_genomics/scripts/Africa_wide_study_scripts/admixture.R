# Load required libraries
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(readr)

# -------------------- FUNCTION DEFINITIONS --------------------

# Function to load and merge ADMIXTURE data with metadata
load_admixture_data <- function(q_file, fam_file, metadata_file) {
  sample_names <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  admix_q <- read.table(q_file, header = FALSE)
  colnames(admix_q) <- paste0("Cluster_", 1:ncol(admix_q))
  admix_data <- cbind(SampleID = sample_names, admix_q)
  
  metadata <- read_tsv(metadata_file, show_col_types = FALSE)
  admix_meta <- left_join(admix_data, metadata, by = "SampleID")
  return(admix_meta)
}

# Function to reshape ADMIXTURE data for plotting
reshape_admixture_long <- function(admix_meta) {
  cluster_cols <- grep("^Cluster_", names(admix_meta), value = TRUE)
  
  admix_long <- melt(
    admix_meta[, c("SampleID", "Region", cluster_cols)],
    id.vars = c("SampleID", "Region"),
    variable.name = "Cluster",
    value.name = "Ancestry"
  )
  
  admix_long$Ancestry <- as.numeric(admix_long$Ancestry)
  admix_long <- admix_long %>%
    arrange(Region) %>%
    mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
  
  return(admix_long)
}

# Function to generate the country color strip
generate_country_strip <- function(admix_meta, admix_long) {
  country_data <- admix_meta %>%
    select(SampleID, Country, Region) %>%
    distinct() %>%
    mutate(SampleID = factor(SampleID, levels = unique(admix_long$SampleID))) %>%
    arrange(SampleID) %>%
    mutate(Country = factor(Country, levels = unique(Country)))
  
  country_colors <- setNames(
    colorRampPalette(brewer.pal(12, "Set3"))(length(unique(country_data$Country))),
    unique(country_data$Country)
  )
  
  strip_plot <- ggplot(country_data, aes(x = SampleID, y = 1, fill = Country)) +
    geom_tile(color = "black", linewidth = 0.1) +
    facet_wrap(~Region, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = country_colors) + 
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.4, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    guides(fill = guide_legend(title = "Country", nrow = 2, byrow = TRUE))
  
  return(strip_plot)
}

# Function to generate the ADMIXTURE bar plot
generate_admixture_plot <- function(admix_long) {
  ggplot(admix_long, aes(x = SampleID, y = Ancestry, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(~Region, scales = "free_x", nrow = 1) +
    scale_y_continuous(labels = percent_format()) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_blank()
    ) +
    labs(
      x = "Samples coloured by country",
      y = "Ancestry Proportion",
      title = paste("ADMIXTURE Plot (K =", length(unique(admix_long$Cluster)), ") Grouped by Region")
    ) +
    guides(fill = "none")
}

# -------------------- EXECUTION --------------------

# File paths
q_file <- "path/to/admixture_ready.3.Q"
fam_file <- "path/to/admixture_pruned_variants_filtered.fam"
metadata_file <- "path/to/metadata.tsv"

# Load and prepare data
admix_meta <- load_admixture_data(q_file, fam_file, metadata_file)
admix_long <- reshape_admixture_long(admix_meta)

# Generate plots
admix_plot <- generate_admixture_plot(admix_long)
country_strip <- generate_country_strip(admix_meta, admix_long)

# Combine plots
final_plot <- admix_plot / country_strip + plot_layout(heights = c(6, 0.3))
print(final_plot)
