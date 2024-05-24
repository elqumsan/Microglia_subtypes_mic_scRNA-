library(tidyverse)
library(cowplot)
library(Seurat)
#library(SeuratData)
library(SeuratWrappers)
#library(Azimuth)
library(patchwork)
library(extrafont)
library(stringr)
library('tidyr')
library("enrichCellMarkers")
source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")

Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
WT_Project <- "WT_data_Microglia"
AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
AZT_Project <-  "AZT_data_Microglia"

WTdata  <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
WT_object <-  CreateSeuratObject(counts = WTdata, project = "Microglia subtypes_WT", min.cells = 3, min.features = 300)
AZT_object <- CreateSeuratObject(counts = AZTdata, project = "Microglia subtypes_AZT", min.cells = 3, min.features = 300)

integrated_object <-  merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), merge.data= TRUE)

integrated_object[["Strain"]] <- factor(integrated_object@meta.data$orig.ident, levels = strain)
integrated_object$Strain <- str_replace(integrated_object$Strain,pattern = "Microglia subtypes_WT", replacement = "WT" )
integrated_object$Strain <- str_replace(integrated_object$Strain, pattern ="Microglia subtypes_AZT", replacement = "AZT" )

integrated_object[["percent.mt"]] <-PercentageFeatureSet(integrated_object, pattern = "^MT-")

integrated_object <- subset(integrated_object, subset = nFeature_RNA > 600 & percent.mt < 8)




#GetAssayData(integrated_object@assays$RNA)

#integrated_object[["joind"]] <- Jointlyers(integrated_object, assay = "RNA")

####### ribosomal gene

ribo.genes <- grep(pattern = "^Rp[s1][[:digit:]]", x = rownames(integrated_object@assays$RNA), value = TRUE)
integrated_object$percent.ribo <- PercentageFeatureSet(integrated_object,features = ribo.genes)



####### All myeloid cells 
path= paste("../Microglia_subtypes_mic_scRNA-/findings/01a_QC_strains/", "/",sep = "")

integrated_object <- integrated_object %>% NormalizeData() 
integrated.strain <- subset(integrated_object ) %>% NormalizeData()

integrated.strain <- AddModuleScore(object = integrated.strain@assays$RNA, features = ribo.genes)

integrated.strain<- integrated.strain %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("batch","ribo.genes", "percent.mt", "nFeature_RNA")) %>%
  RunPCA()
integrated.strain <- JackStraw(integrated.strain, num.replicate = 30, dims = 30)
integrated.strain<-ScoreJackStraw(integrated.strain,dims = 1:30)

ElbowPlot(integrated.strain,ndims = 30) + ggtitle(label = paste("AZT", sep = ""))
ggsave(paste(path,"ElbowPlot", ".png", sep = "" ), width = 7, height = 4, dpi = 150)
print(integrated.strain[["pca"]], dims = 1:30, nfeatures = 30)

# Choose dimensions of PCA after checking ElbowPlot and positive genes in each Clusters
pca_dim <- 11
integrated.strain <- integrated.strain %>%
      RunUMAP(reduction = "pca", dims = 1:pca_dim) %>%
      FindNeighbors(reduction = "pca", dims = 1:pca_dim) %>%
      FindClusters(resolution = 0.5)

res="res05" # specify which rresolution used to cluster cells

# check dimension reduction

DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 0.001) +
  ggtitle(label = "AZT & WT") + coord_fixed()
ggsave(paste(path, "DimPlot1",".png", sep = ""), units = "in", width= 7.3 , height = 7, dpi= 150)

DimPlot(integrated.strain, reduction = "umap" ,label = FALSE, group.by = "Strain", pt.size = 0.001 )
ggsave(paste(path,"_", res, "DimPlot3", ".png", sep = ""), units = "in", width = 7.3, height = 7, dpi = 150)

# DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 0.001, split.by = "seurat_clusters")

QC_plot(integrated.strain@meta.data)
