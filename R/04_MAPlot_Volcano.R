# 05_MAPlot_Volcano.R
#' run_MA_Plots: Generate MA plots for multiple comparisons and save as a grid
#'
#' @param res Named list of data.frames of ranked DESeq2 results (with columns gene_name, log2FoldChange, padj)
#' @param output_file Path to save the MA plots PNG (default: file.path("Results", "MAplot.png"))
#' @return Invisibly returns NULL
run_MA_Plots <- function(res,
                         output_file = file.path("Results", "MAplot.png")) {
  dir.create("Results", showWarnings = FALSE)
  library(ggpubr)
  library(cowplot)

  ma_plot <- function(df) {
    ggmaplot(
      df,
      fc = 2,
      size = 1,
      palette = c("#B31B21", "#1465AC", "darkgray"),
      genenames = as.vector(df$gene_name),
      legend = "bottom",
      top = 20,
      font.label = c("bold", 12),
      label.rectangle = TRUE,
      font.legend = "bold",
      font.main = "bold",
      ggtheme = theme_classic()
    )
  }

  ggsave(filename = output_file,
         plot = ma_plot(res),
         width = 6,
         height = 6,
         units = "in",
         dpi = 300)
  invisible(NULL)
}

#' run_Volcano_Plots: Generate Volcano plots for multiple comparisons and save as a grid
#'
#' @param res Named list of data.frames of ranked DESeq2 results (with columns gene_name, log2FoldChange, padj)
#' @param output_file Path to save the Volcano plots PNG (default: file.path("Results", "Volcanoplot.png"))
#' @return Invisibly returns NULL
run_Volcano_Plots <- function(res,
                              output_file = file.path("Results", "Volcanoplot.png")) {
  dir.create("Results", showWarnings = FALSE)
  library(EnhancedVolcano)

  volcano_plot <- function(df) {
    keyvals <- ifelse(df$log2FoldChange > 1 & df$padj < 0.05, 'red',
               ifelse(df$log2FoldChange < -1 & df$padj < 0.05,'royalblue','grey'))
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == 'red'] <- 'Upregulated'
    names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'
    
    labels <- df %>%
      mutate(score = abs(log2FoldChange) + abs(log10(padj))) %>%
      filter(abs(log2FoldChange) > 1) %>%
      arrange(desc(score)) %>%
      slice_head(n = 15) %>%
      pull(gene_name)

    p <- EnhancedVolcano(
      df,
      lab = df$gene_name,
      selectLab = labels,
      title = "",
      subtitle = "",
      x = 'log2FoldChange',
      y = 'padj',
      xlab = bquote(~Log[2]~"fold change"),
      pCutoff = 0.05,
      FCcutoff = 1,
      colCustom = keyvals,
      pointSize = 2.0,
      labSize = 6.0,
      boxedLabels = TRUE,
      drawConnectors = TRUE,
      widthConnectors = 1.0,
      colConnectors = 'black',
      legendPosition = 'right',
      legendLabSize = 14
    )
    return(p)
  }

  ggsave(filename = output_file,
         plot = volcano_plot(res),
         width = 12,
         height = 12,
         units = "in",
         dpi = 300)
  invisible(NULL)
}
