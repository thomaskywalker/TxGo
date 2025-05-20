# 03_Heatmap.R ------------------------------------------------------------
# Draw a variance-stabilised expression heat-map directly from a DESeq2 object
# (c) 2025 – Shen Lab Pipeline

#' run_Heatmap
#'
#' @description
#' Produce a sample-annotated heat-map from the variance-stabilised matrix
#' (`assay(vst(dds))`).  You may optionally remove a batch effect and/or keep
#' only the *n* most variable genes.  The PNG is written to *Results/*.
#'
#' @param dds          A `DESeq2` object that has already been run through `DESeq()`
#' @param coldata      Sample metadata (**rownames = samples**) containing a
#'                     column named **Group**
#' @param group_colors Named vector mapping Group levels to colours
#' @param batch_col    (string | NULL)  Name of a metadata column to remove as
#'                     batch effect; `NULL` = skip batch correction (default `NULL`)
#' @param top_n        (int | NULL) Keep the *n* most variable genes;
#'                     `NULL` = plot all genes (default `NULL`)
#' @param clustering_distance_rows Distance metric for rows
#'                     (`"correlation"`, `"euclidean"`, …).  Default `"correlation"`
#' @param clustering_distance_cols Distance metric for columns.  Default `"correlation"`
#' @param clustering_method        Agglomeration method (`"ward.D2"`, `"complete"`, …)
#' @param angle_col                Angle for column labels (default `45`)
#' @param scale                    `"row"`, `"column"`, or `"none"` (default `"row"`)
#' @param fontsize                 Base font size (default `10`)
#' @param output_file              PNG path (default `"Results/heatmap.png"`)
#'
#' @return `NULL` (invisible).  The heat-map is saved to `output_file`.
#' @export
#'
run_Heatmap <- function(dds,
                        coldata,
                        group_colors,
                        batch_col                = NULL,
                        top_n                    = NULL,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        clustering_method        = "ward.D2",
                        angle_col                = 45,
                        scale                    = "row",
                        fontsize                 = 10,
                        output_file              = file.path("Results", "heatmap.png")) {

  ## --------------------------------------------------------------------- ##
  ## 0  Ensure output directory exists
  dir.create("Results", showWarnings = FALSE)

  ## --------------------------------------------------------------------- ##
  ## 1  Load required packages quietly
  suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
    library(matrixStats)
    library(pheatmap)
  })

  ## --------------------------------------------------------------------- ##
  ## 2  Variance-stabilising transformation
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)

  ## --------------------------------------------------------------------- ##
  ## 3  Optional batch-effect removal
  if (!is.null(batch_col) && batch_col %in% colnames(coldata)) {
    mat <- removeBatchEffect(mat, batch = coldata[[batch_col]])
  }

  ## --------------------------------------------------------------------- ##
  ## 4  Optional: keep the top n most variable genes
  if (!is.null(top_n)) {
    vars <- rowVars(mat)
    mat  <- mat[
      order(vars, decreasing = TRUE)[seq_len(top_n)],
      ,
      drop = FALSE
    ]
  }

  ## --------------------------------------------------------------------- ##
  ## 5  Build column annotation
  ann <- data.frame(
    Group = factor(coldata$Group, levels = names(group_colors)),
    row.names = rownames(coldata)
  )
  ann_colors <- list(Group = group_colors)

  colnames(ann) <- "Group          "
  names(ann_colors) <- "Group          "
  
  ## --------------------------------------------------------------------- ##
  ## 6  Plot and save
  pheatmap(
    mat,
    annotation_col           = ann,
    annotation_colors        = ann_colors,
    show_rownames            = FALSE,
    show_colnames            = TRUE,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_cols = clustering_distance_cols,
    clustering_method        = clustering_method,
    angle_col                = angle_col,
    scale                    = scale,
    fontsize                 = fontsize,
    filename                 = output_file,
    silent                   = TRUE
  )

  invisible(NULL)
}
