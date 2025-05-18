# 02_PCA.R
#' run_PCA: Perform PCA on DESeq2 output, optionally remove batch effects, and save a biplot in Results/
#'
#' @param dds DESeqDataSet object returned by run_DESeq2()
#' @param coldata Data.frame of sample metadata with rownames matching dds column names; must include a "Group" column and optionally a batch column
#' @param group_names Named list with elements "treatment" and "control" corresponding to levels in coldata$Group
#' @param group_colors Named character vector mapping Group levels to colors
#' @param batch_col Character name of batch column in coldata (default: NULL, no batch correction)
#' @param var_threshold Numeric proportion of variance to remove (default: 0.1)
#' @param output_file Path to save the PCA plot, default in Results/ (default: file.path("Results","PCA.png"))
#' @return A list with elements:
#'   - rld: the variance-stabilized, batch-corrected matrix
#'   - plot: ggplot object of the PCA biplot
run_PCA <- function(dds,
                    coldata,
                    group_names,
                    group_colors,
                    batch_col    = NULL,
                    var_threshold = 0.1,
                    output_file  = file.path("Results","PCA.png")) {
  # Ensure Results directory exists
  dir.create("Results", showWarnings = FALSE)

  # Load required packages
  library(PCAtools)
  library(limma)
  library(ggplot2)

  # 1. Variance-stabilizing transformation
  vsd_mat <- assay(vst(dds, blind = FALSE))

  # 2. Optional batch effect removal
  if (!is.null(batch_col) && batch_col %in% colnames(coldata)) {
    vsd_mat <- limma::removeBatchEffect(vsd_mat,
                                         batch = coldata[[batch_col]])
  }

  # 3. Align sample order
  common <- intersect(colnames(vsd_mat), rownames(coldata))
  vsd_mat <- vsd_mat[, common, drop = FALSE]
  coldata <- coldata[common, , drop = FALSE]

  # 4. Ensure Group factor levels match group_names
  coldata$Group <- factor(coldata$Group,
                          levels = c(group_names$treatment,
                                     group_names$control))

  # 5. PCA analysis
  pca_res <- PCAtools::pca(vsd_mat,
                           metadata    = coldata,
                           removeVar   = var_threshold)

  # 6. Create biplot
#   p <- PCAtools::biplot(pca_res,
#                         x             = "PC1",
#                         y             = "PC2",
#                         colby         = "Group",
#                         colkey        = group_colors,
#                         lab           = NULL,
#                         titleLabSize  = 16,
#                         legendPosition = "right")
  p <- PCAtools::biplot(pca_res,
                        x             = "PC1",
                        y             = "PC2",
                        colby         = "Group",
                        colkey        = group_colors,
                        labSize = 4,
                        encircle = TRUE,
                        encircleFill = TRUE,
                        hline = 0, vline = 0,
                        legendPosition = "right", legendLabSize = 16, legendIconSize = 6.0)
  # 7. Save plot
  ggsave(filename = output_file,
         plot     = p,
         width    = 9,
         height   = 8,
         units    = "in",
         dpi      = 300)

  return(list(rld  = vsd_mat,
              plot = p))
}
