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

## Installation

```bash
# clone the repository
git clone https://github.com/thomaskywalker/TxGo.git
cd TxGo
```

## How to run (Run in VScode or use R project)
1. Open **`TxGo.R`**
2. Run `source("R/Run.R") `
3. Add you data to the `data/` folder
4. Update path to your count data and coldata 
5. Finish the configs
6. Run analysis
7. Wait and see the results

### Hint: You can always check the demo data before using yours