---
title: "Create Seurat object"
output: 
  html_document:
    keep_md: true
---


```r
knitr::opts_chunk$set(eval = FALSE, error = FALSE, eval= FALSE)
```


```r
# Suppress Loading message when build the report 
suppressPackageStartupMessages(
  {
    
library(recount)
library(SummarizedExperiment)
library(S4Vectors)

library(knitr)
library(SummarizedExperiment)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork) 
library(BiocManager)
library(glmGamPoi)
    
  }
)
```

# Loading datasets 


```r
Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
WT_Project <- "WT_data_Microglia"
AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
AZT_Project <-  "AZT_data_Microglia"

WTdata <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
```


# Creat Seurat objects


```r
WT_object <- CreateSeuratObject(counts = WTdata , project = WT_Project , min.cells = 3, min.features = 300)
AZT_object <- CreateSeuratObject(counts = AZTdata , project = AZT_Project , min.cells = 3, min.features = 300)
```

# Merge individual Seurat Objects into Single Seurat Object


```r
merged_objects <- merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), project = "Merged WT_AZT", 
                            merge.data = TRUE )
```
