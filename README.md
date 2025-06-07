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

## Usage

1. Launch the Shiny app directly from GitHub:

```r
library(shiny)
runGitHub(repo = "seurat_visualisation_shiny_app", user = "robban13", ref = "main")
```

2. A browser window will open. Upload an already clustered Seurat object when prompted. The object should have normalisation, clustering and UMAP computed. If your data contains multiple tissues or samples, ensure a metadata column specifies the tissue of each cell.

3. Explore the different tabs (UMAP, Feature Plot, Volcano Plot, Heatmap, Spatial Plot, etc.) to visualise and analyse your data. Use the sidebar to group tissues, rename or merge clusters and download results. The app allows uploads up to 6 GB by default.
