############################################################
# Script 03: SourceTracker-based source–sink analysis
#
# Purpose:
#   This script estimates the contribution of early-stage
#   (Bud) microbiomes as sources to later developmental
#   stage microbiomes (Flower, R3, R5) using SourceTracker v1
#   (Knights et al.).
#
#   Source–sink inference is performed based on integer ASV
#   count tables derived from DADA2 processing.
#
# Major steps:
#   1. Load and validate SourceTracker dependencies
#   2. Import ASV count table and sample metadata
#   3. Define source (Bud) and sink (Flower/R3/R5) samples
#   4. Train SourceTracker model on source samples
#   5. Predict source contributions in sink samples
#   6. Summarize Bud and Unknown contributions by stage
#   7. Generate publication-ready visualizations
#
# Input:
#   - ASV count table (integer counts; rows = samples, cols = ASVs)
#   - Sample metadata with Stage information
#
# Output:
#   - Source contribution table for all sink samples
#   - Bud-only contribution table
#   - Stage-wise mean ± SD summaries
#   - Stacked bar plots (PNG/PDF)
#
# This script should be executed after:
#   Script 01: DADA2-based ASV inference
#   Script 02: Community ecology analysis (microeco)
#
# Author: [Your Name]
# Affiliation: [Your Institution]
# Date: [YYYY-MM-DD]
############################################################


############################
# 1. Load required packages
############################

required_pkgs <- c("vegan", "MCMCpack", "dplyr", "tidyr", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(vegan)
library(MCMCpack)
library(dplyr)
library(tidyr)
library(ggplot2)


############################
# 2. Load SourceTracker
############################

# SourceTracker v1 is not an R package and must be sourced
# from the original GitHub repository (Knights et al.)
#
# https://github.com/danknights/sourcetracker

cat("Select the SourceTracker directory (containing src/SourceTracker.r)\n")
SOURCETRACKER_DIR <- choose.dir(caption = "Select sourcetracker directory")

if (is.na(SOURCETRACKER_DIR)) {
  stop("SourceTracker directory selection was cancelled.")
}

st_file <- file.path(SOURCETRACKER_DIR, "src", "SourceTracker.r")
if (!file.exists(st_file)) {
  stop("src/SourceTracker.r not found in the selected directory.")
}

setwd(SOURCETRACKER_DIR)
source(st_file)


############################
# 3. Define input file paths
############################

# ASV count table (must contain integer counts)
otu_fp  <- "~/Large_Source_tracker/cpm_normalization.csv"

# Sample metadata (must contain Stage information)
meta_fp <- "~/Large_Source_tracker/SAM_info.csv"


############################
# 4. Import and format data
############################

otu_raw  <- read.csv(otu_fp,  check.names = FALSE, stringsAsFactors = FALSE)
meta_raw <- read.csv(meta_fp, check.names = FALSE, stringsAsFactors = FALSE)

# Extract ASV identifiers
feature_ids <- otu_raw$ASVs

# Identify sample columns
sample_cols <- setdiff(
  colnames(otu_raw),
  c("ASVs", "taxonomic_name", "taxa_level")
)

# Convert to sample × feature matrix
otu_mat <- as.matrix(otu_raw[, sample_cols])
rownames(otu_mat) <- feature_ids
otu_mat <- t(otu_mat)

# Align metadata
meta <- meta_raw
rownames(meta) <- meta$name

# Keep only shared samples
common_ids <- intersect(rownames(meta), rownames(otu_mat))
if (length(common_ids) <= 1) {
  stop("Too few shared samples between metadata and ASV table.")
}

otu_mat <- otu_mat[common_ids, , drop = FALSE]
meta    <- meta[common_ids, , drop = FALSE]


############################
# 5. Define source and sink samples
############################

# Bud samples are treated as sources
# Flower, R3, and R5 samples are treated as sinks
meta$SourceSink <- ifelse(meta$Stage == "Bud", "source", "sink")

# SourceTracker uses environment labels for training
meta$Env <- ifelse(meta$Stage == "Bud", "Bud", meta$Stage)

train_ix <- which(meta$SourceSink == "source")
test_ix  <- which(meta$SourceSink == "sink")

envs <- meta$Env


############################
# 6. Validate input matrix
############################

# SourceTracker v1 requires integer count data
if (any(otu_mat < 0)) {
  stop("Negative values detected in ASV table.")
}

if (any(abs(otu_mat - round(otu_mat)) > 1e-9)) {
  stop("ASV table must contain integer counts (not relative abundances).")
}


############################
# 7. Train SourceTracker model
############################

# Use default alpha parameters (recommended for large datasets)
alpha1 <- 0.001
alpha2 <- 0.001

st <- sourcetracker(
  otu_mat[train_ix, , drop = FALSE],
  envs[train_ix]
)


############################
# 8. Predict source contributions in sink samples
############################

results_sink <- predict(
  st,
  otu_mat[test_ix, , drop = FALSE],
  alpha1 = alpha1,
  alpha2 = alpha2
)

# Extract mixture proportions
mix_df <- as.data.frame(results_sink$proportions)
mix_df$SampleID <- rownames(mix_df)
mix_df$Stage <- meta[rownames(mix_df), "Stage"]


############################
# 9. Summarize Bud contribution
############################

if (!("Bud" %in% colnames(mix_df))) {
  stop("Bud column not found in SourceTracker output.")
}

bud_only <- mix_df %>%
  select(SampleID, Stage, Bud) %>%
  arrange(Stage)

bud_summary <- bud_only %>%
  group_by(Stage) %>%
  summarise(
    MeanBud = mean(Bud),
    SDBud   = sd(Bud),
    .groups = "drop"
  )


############################
# 10. Save results
############################

write.csv(
  mix_df,
  "~/Large_Source_tracker/SourceTracker_sink_mixture_full.csv",
  row.names = FALSE
)

write.csv(
  bud_only,
  "~/Large_Source_tracker/SourceTracker_sink_Bud_only.csv",
  row.names = FALSE
)

write.csv(
  bud_summary,
  "~/Large_Source_tracker/SourceTracker_Bud_contribution_by_stage.csv",
  row.names = FALSE
)


############################
# 11. Visualization: stacked bar plots
############################

# Identify Unknown column
mix_cols <- setdiff(colnames(mix_df), c("SampleID", "Stage"))
unknown_col <- mix_cols[grepl("unknown|unassigned", mix_cols, ignore.case = TRUE)][1]

plot_long <- mix_df %>%
  filter(Stage %in% c("Flower", "R3", "R5")) %>%
  pivot_longer(
    cols = c("Bud", unknown_col),
    names_to = "Source",
    values_to = "Proportion"
  )

p <- ggplot(plot_long, aes(x = SampleID, y = Proportion, fill = Source)) +
  geom_col() +
  facet_grid(. ~ Stage, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Samples",
    y = "Proportion",
    title = "SourceTracker-estimated Bud contributions to sink stages"
  )

ggsave("Sink_Bud_Unknown_stacked_barplot.png", p, width = 14, height = 5, dpi = 300)
ggsave("Sink_Bud_Unknown_stacked_barplot.pdf", p, width = 14, height = 5)


############################################################
# End of Script 03
############################################################
