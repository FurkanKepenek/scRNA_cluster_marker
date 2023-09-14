# rds file is available here: https://drive.google.com/file/d/1HfuIqU9NVefos6gZu_SLZUd6GUaqlaW6/view?usp=drive_link
# this file has been gathering from here: https://github.com/FurkanKepenek/singlecell_harmony_R 
# scripts for identifying clusters and markers in single cell RNA-seq data

# install needed packages

library(Seurat)
library(tidyverse)
library(metap)
library(multtest)

# load data to system, you can chose directory for data

ifnb_harmony_data <- readRDS('ifnb_harmony_data.rds')

# get general data information, as understand Seurat is a S4 object and there are sub parts under it 

str(ifnb_harmony_data)

# observe meta data part of Seurat object

View(ifnb_harmony_data@meta.data)

# ---------------------------------Visualization--------------------------------

# firstly clusters will be visualized, reduction parameter should be umap because 
# UMAP algorithm used before 

cluster_plot <- DimPlot(ifnb_harmony_data, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)

# there are two condition for data, control (CTRL) and treated (STIM),
# visualization stim has been choosen

condition_plot <- DimPlot(ifnb_harmony_data, reduction = 'umap', group.by = 'stim')

# compare the plots, first plot differentiates the cells according to condition
# STIM or CTRL, second plot shows cell clusters

condition_plot|cluster_plot

# -----------------------------Finding All Markers------------------------------

# for this workflow there is no need to find all markers but it 
# should be useful in the future

FindAllMarkers(ifnb_harmony_data,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# in this part specific conserved markers will be identified 
# main data (shared with drive) is already normalized and pre processed 
# also, to learn the type of assay for data DefaultAssay() function will be used


DefaultAssay(ifnb_harmony_data)

# to get specific class of cluster FindConservedMarkers() function used
# in the first argument data defined after then we decide class of cluster to 
# define markers for. Finally we chose the group as a "stim"
# basically there will be a comparison between cluster 3 and other clusters

cluster_markers_3 <- FindConservedMarkers(ifnb_harmony_data,
                     ident.1 = 3,
                     grouping.var = 'stim')

# observe the result generally

head(cluster_markers_3)

# Visualization of top feature, first gene from "cluster_markers_3" data frame

FeaturePlot(ifnb_harmony_data, features = c('FCGR3A'), min.cutoff = 'q10')


# renaming cluster 3 name, it is known before please have a look at harmony 
# integration repository :)
# observe the name of current indents, from 0 to 13

Idents(ifnb_harmony_data)

# select main data, choose the cluster after then give the new name of cluster 

ifnb_harmony_data <- RenameIdents(ifnb_harmony_data, `3` = 'CD16 Mono')

# create dim plot with renamed cluster
 
DimPlot(ifnb_harmony_data, reduction = 'umap', label = T)

# all of the cells already have annotation, when meta data downloaded from 
# "SeuratData" package before it comes with annotation
# it can be seen as a column which "seurat_annotations"

View(ifnb_harmony_data@meta.data)


# change other ident names from annotation

Idents(ifnb_harmony_data) <- ifnb_harmony_data@meta.data$seurat_annotations

# observe the name of current indents, from 0 to 13

Idents(ifnb_harmony_data)

# create dim plot with renamed clusters

DimPlot(ifnb_harmony_data, reduction = 'umap', label = TRUE)


# find markers according to conditions "CTRL" and "STIM"
# add CTRL and STIM label to idents

ifnb_harmony_data$celltype.cnd <- paste0(ifnb_harmony_data$seurat_annotations,'_', ifnb_harmony_data$stim)

# observe labelled data

View(ifnb_harmony_data@meta.data)

# rename ıdents according to labelled cell names
Idents(ifnb_harmony_data) <- ifnb_harmony_data$celltype.cnd

DimPlot(ifnb_harmony_data, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(ifnb_harmony_data, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(cluster_markers_3)


FeaturePlot(ifnb_harmony_data, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')



