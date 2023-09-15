## scRNA_cluster_marker

```markdown
# RDS file is available [here](https://drive.google.com/file/d/1HfuIqU9NVefos6gZu_SLZUd6GUaqlaW6/view?usp=drive_link).
# This file has been gathered from [here](https://github.com/FurkanKepenek/singlecell_harmony_R).
# Scripts for identifying clusters and markers in single-cell RNA-seq data.

## Install Needed Packages

```R
library(Seurat)
library(tidyverse)
library(metap)
library(multtest)
```

## Load Data to System

You can choose the directory for data.

```R
ifnb_harmony_data <- readRDS('ifnb_harmony_data.rds')
```

## Get General Data Information

As I understand, Seurat is an S4 object and there are subparts under it.

```R
str(ifnb_harmony_data)
```

## Observe Meta Data Part of Seurat Object

```R
View(ifnb_harmony_data@meta.data)
```

## Visualization

### Cluster Visualization

First, clusters will be visualized. The reduction parameter should be "umap" because UMAP algorithm was used before.

```R
cluster_plot <- DimPlot(ifnb_harmony_data, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
```

### Condition Visualization

There are two conditions for data: control (CTRL) and treated (STIM). Visualization for STIM has been chosen.

```R
condition_plot <- DimPlot(ifnb_harmony_data, reduction = 'umap', group.by = 'stim')
```

Compare the plots: the first plot differentiates the cells according to the condition (STIM or CTRL), and the second plot shows cell clusters.

```R
condition_plot | cluster_plot
```

## Finding All Markers

For this workflow, there is no need to find all markers, but it should be useful in the future.

```R
FindAllMarkers(ifnb_harmony_data,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')
```

## Specific Conserved Markers

In this part, specific conserved markers will be identified. The main data (shared with the drive) is already normalized and pre-processed. Also, to learn the type of assay for data, the `DefaultAssay()` function will be used.

```R
DefaultAssay(ifnb_harmony_data)
```

To get a specific class of cluster, the `FindConservedMarkers()` function is used. In the first argument, data is defined, and then we decide the class of cluster to define markers for. Finally, we chose the group as "stim." Basically, there will be a comparison between cluster 3 and other clusters.

```R
cluster_markers_3 <- FindConservedMarkers(ifnb_harmony_data,
                     ident.1 = 3,
                     grouping.var = 'stim')
```

Observe the result generally.

```R
head(cluster_markers_3)
```

Visualization of the top feature, the first gene from the "cluster_markers_3" data frame.

```R
FeaturePlot(ifnb_harmony_data, features = c('FCGR3A'), min.cutoff = 'q10')
```

## Renaming Cluster 3 Name

It is known before, please have a look at the harmony integration repository. Observe the name of current idents, from 0 to 13.

```R
Idents(ifnb_harmony_data)
```

Select the main data, choose the cluster, and then give the new name of the cluster.

```R
ifnb_harmony_data <- RenameIdents(ifnb_harmony_data, `3` = 'CD16 Mono')
```

Create a dim plot with renamed cluster.

```R
DimPlot(ifnb_harmony_data, reduction = 'umap', label = TRUE)
```

All of the cells already have annotations. When metadata was downloaded from the "SeuratData" package before, it comes with annotation. It can be seen as a column which is "seurat_annotations."

```R
View(ifnb_harmony_data@meta.data)
```

Change other ident names from annotation.

```R
Idents(ifnb_harmony_data) <- ifnb_harmony_data@meta.data$seurat_annotations
```

Observe the name of current idents, from 0 to 13.

```R
Idents(ifnb_harmony_data)
```

Create a dim plot with renamed clusters.

```R
DimPlot(ifnb_harmony_data, reduction = 'umap', label = TRUE)
```

## Find Markers According to Conditions

Find markers according to conditions "CTRL" and "STIM." Add CTRL and STIM labels to idents.

```R
ifnb_harmony_data$celltype.cnd <- paste0(ifnb_harmony_data$seurat_annotations,'_', ifnb_harmony_data$stim)
```

Observe labeled data.

```R
View(ifnb_harmony_data@meta.data)
```

Rename idents according to labeled cell names.

```R
Idents(ifnb_harmony_data) <- ifnb_harmony_data$celltype.cnd
```

Create a dim plot with renamed clusters.

```R
DimPlot(ifnb_harmony_data, reduction = 'umap', label = TRUE)
```

## Find Markers

```R
b.interferon.response <- FindMarkers(ifnb_harmony_data, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')
```

```R
head(b.interferon.response)
```

## Plotting Conserved Features vs DE Features Between Conditions

```R
head(cluster_markers_3)
```

```R
FeaturePlot(ifnb_harmony_data, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')
```

Please note that you should run this code in an R environment with the required packages installed and the data file ('ifnb_harmony_data.rds') available in your working directory or provide the correct file path.
```

This Markdown document represents the provided R code with explanations and links. You can copy and paste it into your Markdown editor or documentation for reference.
