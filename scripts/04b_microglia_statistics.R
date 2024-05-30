library(tidyverse)
library(cowplot)
library(Seurat)

### load object

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
  
