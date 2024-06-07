library(tidyverse)
library(cowplot)
library(Seurat)

### load object
DefaultAssay(integrated_object) <- "RNA"
DefaultAssay(integrated.strain) <- "RNA"

sum_table <-  mg.strain@meta.data %>% group_by(seurat_clusters) %>% 
          summarise( N=n(), ave_nCount_RNA=median(nCount_RNA), ave_nFeature_RNA=median(nFeature_RNA), ave_percent.mt=median(percent.mt))
prop.table(table(Idents(integrated_object),integrated_object$Strain), margin = 2)


#### Plot both genotypes in all strains (all replicates combined)
# generate meta data, 
integrated.meta <- mg.strain@meta.data %>%
  mutate(strain = factor(Strain, levels = c("AZT","WT")),
         new_clusters= ifelse(seurat_clusters %in% 0:5, "H" ,as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("1", "2","3","4")))%>%
  group_by(strain , new_clusters) %>%
  summarise(N=n())

ggplot(integrated.meta, aes(y=N, s=Genotype, fill= new_clusters )) +
  geom_bar(stat = "identity", position = "fill", color= "red") + 
  labs(y="Fraction", fill = "Clusters") +
  facet_grid(~ Genotupe) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = c("bold", "bold.italic")),
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

############ Box Plot for all microglia
############ generate meta data, for statistical testing and box plot 
integrated.meta.stat <- integrated.strain@meta.data %>%
                mutate( strain=factor(global_var$global$strain, levels = c("WT", "AZT")),
                        new_clusters=factor(seurat_clusters, levels = c("1","2","3","4","5","6", "7"))) %>%
                group_by(strain, new_clusters) %>%
                summarise(Med_nFeature=median(nFeature_RNA),
                          Med_percent_mt= median(percent.mt),
                          
                          N=n()) %>%
  group_by(Percent= N/sum(N)*100)

integratedPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#f0E442", "#0072B2", "#D55E00")  

integrated.meta.stat %>%
  ggplot(aes(y=Percent, x= strain)) +
  geom_boxplot(outlier.size = 0, alph= 0.5) +
  geom_point(aes(color = strain), position = position_jitterdodge(), alpha=0.8) +
  scale_color_manual(values = integratedPalette) +
  theme_bw() +
  facet_grid(new_clusters ~ strain, scales = "free_y")
  
