
library(recount)
library(SummarizedExperiment)
library(S4Vectors)

library(shiny)
library(bslib)
library(ggExtra)
library(purrr)

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
library("BiocManager")
library("glmGamPoi")
library('enrichCellMarkers')

source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")
source("../Microglia_subtypes_mic_scRNA-/scripts/02_QC_starin_split.R")



Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
WT_Project <- "WT_data_Microglia"
AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
AZT_Project <-  "AZT_data_Microglia"

WTdata  <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
WT_object <-  CreateSeuratObject(counts = WTdata, project = "Microglia subtypes_WT", min.cells = 3, min.features = 300)
AZT_object <- CreateSeuratObject(counts = AZTdata, project = "Microglia subtypes_AZT", min.cells = 3, min.features = 300)

integrated_object <-  merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), merge.data= TRUE)
# object.ref <- subset(integrated_object, orig.ident %in% c("Microglia subtypes_WT" , "Microglia subtypes_AZT"))

strain <- c("Microglia subtypes_WT", "Microglia subtypes_AZT")
cols =  c("#888888", "#00AA00")

integrated_object[["Strain"]] <- factor(integrated_object@meta.data$orig.ident, levels = strain)
integrated_object$Strain <- str_replace(integrated_object$Strain,pattern = "Microglia subtypes_WT", replacement = "WT" )
integrated_object$Strain <- str_replace(integrated_object$Strain, pattern ="Microglia subtypes_AZT", replacement = "AZT" )
meta <- integrated_object@meta.data

meta_tidy <- meta %>% 
  select(orig.ident,nCount_RNA, nFeature_RNA, Strain) %>%
  gather(-orig.ident, -Strain ,key= "QC", value="value" )
#VlnPlot(object.ref, features = c("nFeature_RNA","nCount_RNA" ) )

map(~QC_plot(integrated_object@meta.data))

######### ribosomal gene
grep( pattern = "^Rp[s1][[:digit]]", x = rownames(integrated_object@assays$RNA),value = TRUE)
results = CMenrich(gene.list = c('Sall1', 'Hexb', 'Fcrls', 'Gpr43', 'Cx3cr1', 'Tmem119', 'Trem2', 'P2ry12', 'Mertk', 'Pros1','Siglech'), species = 'mouse' )


##################################
what_dims(object_type = WT_object, path = Wtdata.path )

cells <- WhichCells(object.ref)

CellsMeta = object.ref@meta.data
randomnumbers <- runif(25167, 0.0, 1.1)
CellsMeta["Gene_IDs"] <- randomnumbers
head(CellsMeta)
cellsMetaTrim <- subset(CellsMeta, select = c("Gene_IDs"))
object.ref <-AddMetaData(object.ref, cellsMetaTrim)
head(object.ref)
VlnPlot(object.ref, c("nCount_RNA", "Gene_IDs"), ncol = 2)

#object.ref[["RNA"]] <- split( object.ref[["RNA"]] ,f =  object.ref$seurat_clusters )
#################

object.ref <- NormalizeData(object.ref)
object.ref <- FindVariableFeatures(object.ref) 
object.ref <-ScaleData(object.ref)
object.ref <-RunPCA(object.ref)
object.ref<-FindNeighbors(object.ref, dims = 1:30)
object.ref <-FindClusters(object.ref, resolution = 0.05)

object.ref <-RunUMAP(object.ref, dims = 1:30)

DimPlot(object.ref, group.by = c("orig.ident", "seurat_clusters"))

#meta.all <- integrated_object %>% 
#  dplyr::rename(Integrated_object = "integrated_object")  %>%
#  select(ID_prefix, sample_ID, Cust_ID, Exp_batch)
