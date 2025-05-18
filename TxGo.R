#!/usr/bin/env Rscript

# Activate the programme
source("R/Run.R")

# 1. Specify input files
count_file   <- "data/genes.readcount.mRNA.csv"
coldata_file <- "data/coldata.xlsx"

# 2. Define group labels
#    Keys are internal names; values are the actual group labels in metadata
group_names <- list(
  treatment = "Normal",
  control   = "Patient"
)

# 3. Define colors for each internal group key
#    Will be mapped to actual group labels automatically
group_colors <- c(
  treatment = "#1f78b4",
  control   = "#a02c2c"
)
# Map colors to actual group labels
group_colors <- setNames(group_colors, unlist(group_names))

# 4. Specify species: "hsa" = human, "mmu" = mouse
species_code <- "hsa"

# 5. Run all analyses
txgo_results <- TxGo(
  count_file   = count_file,
  coldata_file = coldata_file,
  group_names  = group_names,
  group_colors = group_colors,
  species_code = species_code
)
