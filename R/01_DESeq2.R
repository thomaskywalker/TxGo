# 01_DESeq2.R
#' run_DESeq2: Read input files, perform DESeq2 analysis, annotate DEGs, and export results
#'
#' @param count_file Path to CSV file of raw counts (genes in rows, samples in columns)
#' @param coldata_file Path to Excel file of sample metadata; must include a 'Group' column and optional 'sample' column
#' @param design Formula for DESeq2 design (default: ~ Group)
#' @param species_code Species code: 'hs' for human (default) or 'mm' for mouse; used for Ensembl dataset
#' @return A list containing:
#'   - dds: DESeqDataSet object after running DESeq()
#'   - res: results object with differential expression statistics
run_DESeq2 <- function(counts,
                       coldata,
                       species_code = "hsa", ...) {
  # Create Results directory
  dir.create("Results", showWarnings = FALSE)

  # Load required libraries
  library(readxl)
  library(writexl)
  library(tibble)
  library(DESeq2)
  library(apeglm)
  library(biomaRt)
  library(dplyr)
  counts <- counts[, rownames(coldata)]
  ## 1 建 DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,     # 這裡用的是已經帶有 col_data$Group 的 data.frame
    design    = ~ Group
  )

  ## 2 預過濾 + DESeq()
  keep <- rowSums(counts(dds)) >= 10
  dds  <- dds[keep, ]
  dds  <- DESeq2::DESeq(dds)

  res <- results(dds)
  # 7. Prepare results table
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "ensembl_id")
  # Remove version suffix from Ensembl IDs
  res_df$ensembl_id <- sub("\\..*", "", res_df$ensembl_id)
  # Sort by log2 fold change descending
  res_df <- res_df %>% arrange(desc(log2FoldChange))

  # 8. Annotate with gene symbol and Entrez ID
  if (species_code == "hsa") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    orgdb <- org.Hs.eg.db
  } else {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    orgdb <- org.Mm.eg.db
  }
  res_df <- res_df %>% 
  mutate(entrez_id = mapIds(orgdb, keys = ensembl_id, column = "ENTREZID", keytype = "ENSEMBL") %>%
                   ifelse(is.na(.) | duplicated(.), as.character(ensembl_id), .), .before = 1) %>%
  mutate(gene_name = mapIds(orgdb, keys = ensembl_id, column = "SYMBOL", keytype = "ENSEMBL") %>%
                   ifelse(is.na(.) | duplicated(.), as.character(ensembl_id), .), .before = 1)

  # 9. Export annotated DEGs
  write_xlsx(res_df,
             path = file.path("Results", "DEG_annotated.xlsx"))

  # Return objects
  return(list(dds = dds, res = res_df))
}
