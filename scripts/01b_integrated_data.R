#### Integration Microglia

library(tidyverse)
library(cowplot)
library(Seurat)

source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")
meta <- readRDS("../Microglia_subtypes_mic_scRNA-/data/Meta_Markers.rds")
