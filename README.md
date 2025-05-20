# TxGo ▸  one-click bulk-RNA-seq analysis

TxGo integrated state-of-the-art transcriptome analysis tools into a wrapper function, and is easy to operate. Perfect for bench scientist as a starter kit.

## What it does

1. Differential gene expression (DGE) analysis
2. PCA & Heatmap
3. Volcano & MA plot
4. GSEA (GO, KEGG, Hallmark) analysis and visualisation

## Install (Terminal)

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

## Input files (put in `data/`)
1. **`genes.readcount.mRNA.csv`** rows = Ensembl IDs • columns = sample IDs

|               | 1   | 2   | 3   | 4   | 5   | 6   |
|---------------|-----|-----|-----|-----|-----|-----|
| ENSG00000187961 | 405 | 202 | 355 | 333 | 461 | 241 |
| ENSG00000187583 |  89 | 112 |  59 | 170 | 200 | 162 |
| ENSG00000187642 |  10 |   0 |   0 |   7 |   6 |   2 |
| ENSG00000188290 | 470 | 114 | 170 | 256 |  59 | 116 |
| ENSG00000187634 | 327 | 264 | 159 | 154 | 215 |  73 |
| ... | ... | ... | ... | ... | ... | ... |

2. **`coldata.xlsx`** | Contain at least two columns `ID` (sample ID) · `Group` (Treat/Control) | optional `batch` |

| ID | Group   |
|----|---------|
| 1  | Normal  |
| 2  | Normal  |
| 3  | Normal  |
| 4  | Patient |
| 5  | Patient |
| 6  | Patient |

## Hint: You can always check the demo data before using yours