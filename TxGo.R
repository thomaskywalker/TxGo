#!/usr/bin/env Rscript

# 0. Activate the programme
source("R/Run.R")

# 1. Specify input files
count_file   <- "data/genes.readcount.mRNA.csv"
coldata_file <- "data/coldata.xlsx"

# 2. Define group labels
group_names <- list(
  control   = "Normal",
  treatment = "Patient"
)

# 3. Define colors for each internal group key
group_colors <- c(
  control   = "#2c41a0",
  treatment = "#b41f1f"
)

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
