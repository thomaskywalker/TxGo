#!/usr/bin/env Rscript
# Run.R — install/load deps quietly, then define TxGo()

# 0. 顯示進度開頭
message("Preparing environment\nIt may take a while...")

# 1. 開啟一個暫存檔，把 stdout 和 message 都導到裡面
tmp <- tempfile()
out <- file(tmp, open = "wt")
sink(out)
sink(out, type = "message")


# 0. BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 1. Package lists
cran_pkgs <- c(
  # tidy data / IO
  "tidyverse",   # includes ggplot2, dplyr, etc.
  "readxl",
  "writexl",
  # visualisation helpers
  "ggpubr",
  "cowplot",
  "pheatmap",
  "matrixStats",
  "VennDiagram",
  # gene-set helper
  "msigdbr"
)

bioc_pkgs <- c(
  # core differential-expression + downstream
  "DESeq2",
  "limma",
  # annotation / ID mapping
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  # enrichment & visualisation
  "clusterProfiler",
  "DOSE",
  "enrichplot",
  # high-level plots
  "EnhancedVolcano",
  "PCAtools"
)

# 2. Install + load CRAN packages ----------------------------------------------
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg,
                     repos      = "https://cloud.r-project.org",
                     dependencies = TRUE,
                     quiet      = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE)
  )
}

# 3. Install + load Bioconductor packages --------------------------------------
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE)
  )
}

# 6. 把輸出導回 console
sink(type = "message")
sink()
close(out)
unlink(tmp)

r_scripts <- list.files("R", pattern="^\\d+_.*\\.R$", full.names=TRUE)
invisible(lapply(r_scripts, source))

# ─── 2. 定義主流程函式 TxGo() ───────────────────────────────────
TxGo <- function(count_file,
                 coldata_file,
                 group_names,
                 group_colors,
                 species_code){
  # (a) 讀入 count & metadata
  counts  <- read.csv(count_file, row.names=1, check.names=FALSE)
  coldata <- readxl::read_excel(coldata_file) %>%
    as.data.frame()

  if (!all(c("ID","Group") %in% colnames(coldata))) {
    stop("coldata.xlsx 必須包含 'ID' 與 'Group' 兩個欄位")
  }

  rownames(coldata) <- coldata$ID

  coldata$Group <- factor(
    coldata$Group,
    levels = unlist(group_names)
  )
  
  # (b) 1. DESeq2
  dds_out <- run_DESeq2(
    counts = counts,
    coldata     = coldata,
    species_code = species_code
  )
  
  # (c) 2. PCA
  pca_out <- run_PCA(
    dds          = dds_out$dds,
    coldata      = coldata,
    group_names  = group_names,
    group_colors = group_colors
  )
  
  # (d) 3. Heatmap
  run_Heatmap(
    dds          = dds_out$dds,
    coldata      = coldata,
    group_colors = group_colors,
    batch_col    = NULL,  # or NULL if no batch column
    top_n        = NULL      # plot all  genes
  )
  
  # (f) 5. MA & Volcano
  run_MA_Plots(dds_out$res)
  run_Volcano_Plots(dds_out$res)
  
  # (g) 6. GSEA
  gsea_list <- run_GSEA(dds_out$res, species_code = species_code)
  
  # (h) 7. GSEA plot
  run_GSEAplots(res_list = gsea_list)
  
  
  message("=== All analyses are done === \n=== Happy Bioinformatics ===")
}


message("Done!")