



what_dims <- function(object_type , path,  .){
  
  object_type <-AddModuleScore(object = object_type , features = genes, ctrl = 100, name = 'ribo_Features' )
  
  return(object_type)
}


# QC: Violin plot of nfeatures_RNA, percentmt, for ecah cluster

QC_plot <- function(data, y ){
  
  p <- ggplot(data, aes_string(x= "seurat_clusters" , y=y , colors= "seurat_clusters")) + 
    geom_violin()+
    geom_boxplot(width= 0.07, outlier.shape = NA, color = "black", alpha= 0.7) +
    theme_bw()+
    theme(legend.position = "none", axis.title.x =  element_blank())
  return(p)
  
}
  