
library(tidyverse)
library(cowplot)
library(Seurat)
library(patchwork)


Idents(integrated.strain) <- "seurat_clusters"
plot_title = "mic_cluster"
i=17 #(PCA dim)
j=0.6 #(resolution)

##### UMAP Plot

DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 0.001, label.size = 5) +
  coord_fixed() +
  theme(axis.title = element_blank(), legend.position = "none")
ggsave(paste(global_var$global$path_microglia_clustering, "/" ,plot_title, "_", "umap_", i, "_res_", j, "_Dimplot_Strain_small", ".png", sep=""), units = "in", 
       width = 4 , height = 4, dpi = 300 )


##### check marker genes of each cluster, list top 10 marker genes 

integrated_markers %>% group_by(cluster) %>% top_n(-10, wt = p_val_adj)


####### UMAP plot strain  splitted 

Idents(integrated.strain) <- "final_clusters"
p <- DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 1E-10, split.by = "strain", ncol = 2, label.size = 3 , repel = TRUE ) +
  coord_fixed() +
  theme()

ggsave(paste(global_var$global$path_microglia_clustering, "/", "Dimplot_strain_no_text.png", sep = "" ), p , dpi = 300)



######### Dot Plot of microglia marker genes for each clusters 

Idents(integrated.strain) <- "final_clusters"
 genes <- c("Cst7",  "Apoe", "Cx3cr1", "Tmem119", "Ifit3", "Ifitm3", "Irf7", "Hexb", "Cd81", "Cst3", "Rplp1", "Rps21", "Rps24" 
           , "C3ar1", "Stmn1", "Top2a", "Birc5" ) 
 
 file_name <- paste(global_var$global$path_microglia_clustering, "/Dotplot_all_gene.png", sep = "")
 DotPlot(integrated.strain, features = genes) + RotatedAxis() +
   theme(axis.title = element_blank()) +
   scale_y_discrete( labels = function(x) str_wrap(x, width = 20) ) +
   coord_flip() +
   theme(axis.title = element_text( face = "bold.italic"))

 ggsave( file_name, units = "in", width = 5.5 , height = 5.5 , dpi = 300) 
 