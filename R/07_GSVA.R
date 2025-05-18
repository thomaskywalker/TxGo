# 08_GSVA.R
#' run_GSVA: Perform GSVA for Hallmark pathways and save results
#'
#' @param dds DESeqDataSet object returned by run_DESeq2()
#' @param coldata Data.frame of sample metadata with rownames matching dds column names; must include 'Group' and optionally 'batch'
#' @param group_colors Named vector mapping Group levels to colors (for heatmap)
#' @param species_code Species code: 'hs' for human (default) or 'mm' for mouse
#' @return A matrix of GSVA scores (pathways x samples)
run_GSVA <- function(dds,
                     coldata,
                     group_colors,
                     species_code = "hs") {
  # Ensure directories exist
  dir.create("Results", showWarnings = FALSE)
  dir.create(file.path("Results", "GSVA"), showWarnings = FALSE)

  # Load required libraries
  library(GSVA)
  library(msigdbr)
  library(limma)
  library(pheatmap)
  library(matrixStats)
  library(dplyr)
  library(AnnotationDbi)
  if (species_code == "mmu") library(org.Mm.eg.db) else library(org.Hs.eg.db)

  # 1. VST transformation and optional batch correction
  vsd_mat <- assay(vst(dds, blind = FALSE))
  if ("batch" %in% colnames(coldata)) {
    vsd_mat <- limma::removeBatchEffect(vsd_mat, batch = coldata$batch)
  }

  # 2. Map Ensembl IDs to Entrez IDs
  ens_ids <- sub("\\..*", "", rownames(vsd_mat))
  OrgDb <- if (species_code == "mmu") org.Mm.eg.db::org.Mm.eg.db else org.Hs.eg.db::org.Hs.eg.db
  entrez_ids <- mapIds(OrgDb,
                       keys    = ens_ids,
                       column  = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  # Handle NAs or duplicates by falling back to Ensembl ID
  valid <- !is.na(entrez_ids) & !duplicated(entrez_ids)
  rownames(vsd_mat) <- ifelse(valid, entrez_ids, ens_ids)

  # 3. Load Hallmark gene sets
  species_name <- if (species_code == "mm") "Mus musculus" else "Homo sapiens"
  hallmark_df <- msigdbr(species = species_name, category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    as.data.frame()
  hallmark_list <- split(hallmark_df$entrez_gene, hallmark_df$gs_name)

  # 4. Perform GSVA
  gsva_scores <- gsva(as.matrix(vsd_mat), hallmark_list,
                      kcdf = "Gaussian", min.sz = 15, max.sz = 500,
                      mx.diff = TRUE)

  # 5. Save GSVA scores to Excel
  df_scores <- as.data.frame(gsva_scores) %>%
    rownames_to_column(var = "Pathway")
  writexl::write_xlsx(df_scores,
                      path = file.path("Results/GSVA", "GSVA_Hallmark.xlsx"))

  # 6. Heatmap of all Hallmark pathways
  # Prepare column annotation
  ann <- data.frame(Group = factor(coldata$Group, levels = names(group_colors)))
  rownames(ann) <- rownames(coldata)
  ann_colors <- list(Group = group_colors)

  pheatmap(gsva_scores,
           annotation_col    = ann,
           annotation_colors = ann_colors,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           filename          = file.path("Results/GSVA", "GSVA_Hallmark_heatmap.png"),
           width             = 7,
           height            = 10,
           units             = "in",
           silent            = TRUE)

  return(gsva_scores)
}
