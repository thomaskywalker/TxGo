# 07_GSEA_plot.R  --------------------------------------------------------------
# Draw dot- and ridge-plots for a *single* comparison’s GSEA results

#' run_GSEAplots
#'
#' @param res_list   A list like
#'                   list(GO_BP = gseaResult, KEGG = gseaResult, Hallmark = gseaResult)
#' @param comp_name  A label (e.g. "PN") used in plot titles / file names
#' @param categories Categories to plot (must match names in res_list)
#' @param output_dir Folder to save PNGs (default "Results/GSEA_plot")
#'
#' @return `NULL` (plots are written to disk)
#' @export
run_GSEAplots <- function(res_list,
                          comp_name  = "",
                          categories = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Hallmark"),
                          output_dir = file.path("Results", "GSEA_plot")) {

  ## 0  create folder
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## 1  libs
  suppressPackageStartupMessages({
    library(enrichplot)
    library(ggplot2)
  })

  ## 2  iterate over requested categories
  for (cat in intersect(categories, names(res_list))) {

    gr <- res_list[[cat]]
    if (is.null(gr@result$ID) || nrow(gr@result) == 0) next   # skip empty

    # build title / prefix
    prefix <- if (nzchar(comp_name)) paste0(comp_name, "_") else ""
    title  <- if (nzchar(comp_name)) comp_name else "GSEA"

    ## dot plot ----------------------------------------------------------
    dot <- dotplot(gr, orderBy = "NES", showCategory = 5, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste(title, cat)) +
      theme(plot.title = element_text(size = 14))

    ggsave(
      filename = file.path(output_dir, paste0(prefix, cat, "_dot.png")),
      plot     = dot,
      width    = 6, height = 6, dpi = 300
    )

    ## ridge plot --------------------------------------------------------
    ridge <- ridgeplot(gr, orderBy = "NES", showCategory = 10) +
      ggtitle(paste(title, cat)) +
      theme(plot.title = element_text(size = 14))

    ggsave(
      filename = file.path(output_dir, paste0(prefix, cat, "_ridge.png")),
      plot     = ridge,
      width    = 6, height = 6, dpi = 300
    )
  }

  invisible(NULL)
}
