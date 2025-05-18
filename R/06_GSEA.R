# 06_GSEA.R
#' run_GSEA: Perform GSEA analyses (GO, KEGG, Hallmark) on ranked DE results and save outputs
#'
#' @param res Data.frame of ranked DE results; must include columns 'ensembl_id' and 'log2FoldChange'
#' @param fname  Character string to label output files (e.g., comparison name)
#' @param species_code Species code: 'hsa' for human (default) or 'mm' for mouse
#' @return Invisibly returns a list of GSEA results objects
run_GSEA <- function(res,
                     fname,
                     species_code = "hsa") {
  # Ensure output directory
  dir.create("Results/GSEA", recursive = TRUE, showWarnings = FALSE)

  # Load required packages
  library(clusterProfiler)
  library(msigdbr)
  library(DOSE)
  library(writexl)
  library(AnnotationDbi)

  # Determine OrgDb and KEGG species
  if (species_code == "mmu") {
    OrgDb       <- org.Mm.eg.db::org.Mm.eg.db
    kegg_species <- "mmu"
  } else {
    OrgDb       <- org.Hs.eg.db::org.Hs.eg.db
    kegg_species <- "hsa"
  }

  # Replace NA or duplicates with original
  geneList <- res$log2FoldChange
  names(geneList) <- res$entrez_id
  geneList <- sort(geneList, decreasing = TRUE)

  # 2. GSEA for GO: BP, CC, MF
  gseGO_BP <- clusterProfiler::gseGO(geneList, OrgDb = OrgDb, ont = "BP", pvalueCutoff = 1, verbose = FALSE, eps = 0, nPermSimple = 10000)
  gseGO_CC <- clusterProfiler::gseGO(geneList, OrgDb = OrgDb, ont = "CC", pvalueCutoff = 1, verbose = FALSE, eps = 0, nPermSimple = 10000)
  gseGO_MF <- clusterProfiler::gseGO(geneList, OrgDb = OrgDb, ont = "MF", pvalueCutoff = 1, verbose = FALSE, eps = 0, nPermSimple = 10000)
  # Make readable
  gseGO_BP <- DOSE::setReadable(gseGO_BP, OrgDb = OrgDb, keyType = "ENTREZID")
  gseGO_CC <- DOSE::setReadable(gseGO_CC, OrgDb = OrgDb, keyType = "ENTREZID")
  gseGO_MF <- DOSE::setReadable(gseGO_MF, OrgDb = OrgDb, keyType = "ENTREZID")

  # 3. GSEA for KEGG
  gseKEGG <- clusterProfiler::gseKEGG(geneList, organism = kegg_species, pvalueCutoff = 1, verbose = FALSE, eps = 0, nPermSimple = 10000)
  gseKEGG <- DOSE::setReadable(gseKEGG, OrgDb = OrgDb, keyType = "ENTREZID")


  # 4. Hallmark gene sets
  hallmark <- msigdbr(species = if (species_code == "hsa") "Homo sapiens" else "Mus musculus",
                       category = "H", subcategory = NULL) %>%
    dplyr::select(gs_name, entrez_gene) %>%
    as.data.frame()
  gseH <- clusterProfiler::GSEA(geneList, TERM2GENE = hallmark, pvalueCutoff = 1, verbose = FALSE)
  gseH <- DOSE::setReadable(gseH, OrgDb = OrgDb, keyType = "ENTREZID")

  # 5. Save RDS of results
  saveRDS(list(
    GO_BP     = gseGO_BP,
    GO_CC     = gseGO_CC,
    GO_MF     = gseGO_MF,
    KEGG      = gseKEGG,
    Hallmark  = gseH
  ), file = file.path("Results/GSEA", paste0("GSEA_results.rds")))

  # 6. Export to Excel
  res_list <- list(
    BP       = as.data.frame(gseGO_BP),
    CC       = as.data.frame(gseGO_CC),
    MF       = as.data.frame(gseGO_MF),
    KEGG     = as.data.frame(gseKEGG),
    Hallmark = as.data.frame(gseH)
  )
  write_xlsx(res_list,
             path = file.path("Results/GSEA", paste0("GSEA_results.xlsx")))

  return(list(
    GO_BP     = gseGO_BP,
    GO_CC     = gseGO_CC,
    GO_MF     = gseGO_MF,
    KEGG      = gseKEGG,
    Hallmark  = gseH
  ))
}
