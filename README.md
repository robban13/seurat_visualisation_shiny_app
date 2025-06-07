# Seurat Visualisation Shiny App

This repository contains a Shiny application for interactive visualisation and analysis of [Seurat](https://satijalab.org/seurat/) objects.

## Features

- Upload pre-processed Seurat objects (`.rds`).
- Visualise data with UMAP, Feature plots, Volcano plots, heatmaps, violin and dot plots.
- Perform spatial niche and neighbourhood analysis.
- Manage clusters: sub-cluster, rename, merge or group by tissue.
- Download generated plots and tables or save subsetted objects.

## Requirements

- R (>= 4.0 recommended) and an internet connection for installing packages.
- The following R packages are required:
  `shiny`, `Seurat`, `ggplot2`, `EnhancedVolcano`, `dplyr`, `epitools`,
  `pwr`, `DT`, `MAST`, `DESeq2`, `tibble`, `shinyjs`, `clusterProfiler`,
  `slingshot`, `scales`, `viridis`, `RColorBrewer`, `metap`, `RANN`,
  `reshape2`.

To install them:

```r
install.packages(c("shiny", "Seurat", "ggplot2", "EnhancedVolcano", "dplyr",
                   "epitools", "pwr", "DT", "tibble", "shinyjs", "scales",
                   "viridis", "RColorBrewer", "metap", "RANN", "reshape2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("MAST", "DESeq2", "clusterProfiler", "slingshot"))
```

## Usage

1. Clone the repository:

```bash
git clone <repo-url>
cd seurat_visualisation_shiny_app
```

2. Launch the Shiny app directly from GitHub:

```r
library(shiny)
runGitHub(repo = "seurat_visualisation_shiny_app", user = "robban13", ref = "main")
```

   Alternatively, run the app from a cloned copy:

```r
shiny::runApp("app.R")
```

   or from the command line:

```bash
Rscript app.R
```

3. A browser window will open. Upload an already clustered Seurat object when prompted. The object should have normalisation, clustering and UMAP computed. If your data contains multiple tissues or samples, ensure a metadata column specifies the tissue of each cell.

4. Explore the different tabs (UMAP, Feature Plot, Volcano Plot, Heatmap, Spatial Plot, etc.) to visualise and analyse your data. Use the sidebar to group tissues, rename or merge clusters and download results. The app allows uploads up to 6 GB by default.

## License

This project is licensed under the terms of the GNU General Public License v3.0. See [LICENSE](LICENSE) for full details.
