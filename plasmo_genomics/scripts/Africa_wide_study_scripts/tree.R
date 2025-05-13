# ===================================================
#  Phylogenetic Plot with ggtree for Africa wide data
# ===================================================

# ---- Load Required Libraries ----
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(ggtreeExtra)

# ---- Custom Color Palettes ----
country_colors <- c(
  "Angola" = "#1b9e77", "Cameroon" = "#d95f02", "CAR" = "#7570b3", "Congo" = "#e7298a",
  "Equitorial Guinea" = "#66a61e", "Gabon" = "#e6ab02", "Ghana" = "#a6761d", "Ivory Coast" = "#666666",
  "Kenya" = "#1f78b4", "Liberia" = "#b2df8a", "Malawi" = "#fb9a99", "Mali" = "#fdbf6f",
  "Nigeria" = "#cab2d6", "Nigeria Present" = "#6a3d9a", "Sierra Leone" = "#CDC8B1",
  "Sudan" = "#b15928", "Tanzania" = "#8dd3c7"
)

region_colors <- c(
  "Central_Africa" = "#e41a1c", "West_Africa" = "#377eb8",
  "East_Africa" = "#4daf4a", "North_Africa" = "#e6ab02"
)

# ---- Helper Functions ----

# Format tip labels to "SampleID (Country)"
format_labels <- function(metadata) {
  metadata %>%
    mutate(Label = paste0(SampleID, " (", Country, ")"))
}

# Update tree tip labels based on metadata
update_tree_labels <- function(tree, metadata) {
  tree$tip.label <- sapply(tree$tip.label, function(x) {
    if (x %in% metadata$SampleID) {
      paste0(x, " (", metadata$Country[metadata$SampleID == x], ")")
    } else {
      x
    }
  })
  return(tree)
}

# Prepares region annotation data for geom_fruit
prepare_region_annots <- function(metadata) {
  metadata %>%
    mutate(Label = paste0(SampleID, " (", Country, ")"),
           label = ifelse(
             str_detect(Label, "MD[0-9]+"),
             paste0(str_extract(Label, "MD[0-9]+"), " (", Country, ")"),
             Label
           )) %>%
    select(label, Region_Track = Region)
}

# Generate the base tree plot
build_tree_base <- function(tree, layout = "circular", branch_length = "branch.length") {
  ggtree(tree, layout = layout, branch.length = branch_length, size = 0.6, color = "grey30")
}

# Main plotting function
plot_tree <- function(tree_file, metadata_file,
                      layout = "circular",
                      add_region_ring = TRUE,
                      add_tip_labels = TRUE,
                      title = "Phylogenetic Tree") {

  # Load data
  tree <- read.tree(tree_file)
  metadata <- read_tsv(metadata_file, show_col_types = FALSE)
  metadata <- format_labels(metadata)
  tree <- update_tree_labels(tree, metadata)

  # Create base tree
  p <- build_tree_base(tree, layout = layout)

  # Join metadata to tree data
  p$data <- p$data %>%
    left_join(metadata, by = c("label" = "Label")) %>%
    mutate(label = ifelse(
      str_detect(label, "MD[0-9]+"),
      paste0(str_extract(label, "MD[0-9]+"), " (", str_extract(label, "(?<=\\()[^)]+"), ")"),
      label
    ))

  # Region annotation ring data
  if (add_region_ring) {
    region_annots <- prepare_region_annots(metadata)
    p <- p + geom_fruit(
      data = region_annots,
      geom = geom_tile,
      mapping = aes(y = label, fill = Region_Track),
      width = 0.05,
      offset = ifelse(layout == "circular", 0.15, 0),
      color = NA
    )
  }

  # Add tip labels and country color dummy points
  if (add_tip_labels) {
    p <- p + geom_tiplab(
      aes(label = label, color = Country),
      size = 2.3, offset = 0.05, show.legend = FALSE
    )
  }

  # Add dummy points for legend only
  p <- p + geom_point(
    data = p$data %>% filter(!is.na(Country)),
    aes(x = x, y = y, color = Country),
    shape = 16, size = 2, show.legend = TRUE
  )

  # Final styling
  p + scale_color_manual(values = country_colors, name = "Country") +
    scale_fill_manual(values = region_colors, name = "Region") +
    theme_tree2() +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = title)
}

# ===============================
#        Example Usage
# ===============================

# Change these paths to your own files
tree_path <- "tree.treefile"
metadata_path <- "metadata.tsv"

# Call the plotting function
p <- plot_tree(
  tree_file = tree_path,
  metadata_file = metadata_path,
  layout = "circular",           # Options: "circular", "rectangular", etc.
  add_region_ring = TRUE,        # TRUE to add region annotation ring
  add_tip_labels = TRUE,         # TRUE to show tip labels
  title = "Phylogenetic Tree by Country and Region"
)

# Show the plot
print(p)
