# TxGo ▸  one-click bulk-RNA-seq analysis

**What it does**

1. Differential expression
2. PCA & Heatmap
3. Volcano & MA plot
4. GSEA (GO, KEGG, Hallmark)

## Input files (put in `data/`)
1. **`genes.readcount.mRNA.csv`** | *csv* • rows = Ensembl IDs • columns = sample IDs | raw integer counts |
2. **`coldata.xlsx`** | *xlsx* with at least two columns <br>`ID` (sample ID) · `group` (Treat/Control) | optional `batch` column |
---

## How to run (RStudio or VScode)

1. Open **`TxGo.R`**  
2. Finish the configs
3. Run analysis
4. Wait and see the results

## Structure
project/
├── data/
│ ├── genes.readcount.mRNA.csv # raw counts (genes × samples)
│ └── coldata.xlsx # sample metadata: ID | group | (batch)
├── R/ # modular analysis functions
│ ├── 01_DESeq2.R 05_MAPlot_Volcano.R
│ ├── 02_PCA.R 06_GSEA.R
│ ├── 03_Heatmap.R 07_GSEA_plot.R
│ └── 04_VennDiagram.R 08_GSVA.R
├── Run.R # environment setup + TxGo() definition
└── TxGo.R # user config → Source = full analysis