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
  mutate(strain = factor(strain, levels = c("AZT","WT")),
         new_clusters= ifelse(seurat_clusters %in% 0:0, "H" ,as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("1", "2","3","4", "5", "6", "7"))) %>%
  group_by(strain , new_clusters) %>%
  arrange(strain) %>%
  summarise(N=n())

p <- ggplot(integrated.meta, aes(y=N, x=strain , fill= new_clusters )) +
      geom_bar(stat = "identity", position = "fill", color= "red") + 
      labs(y="Fraction", fill = "Clusters") +
      facet_grid(strain~ . ) +
      coord_flip() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x =  element_blank(),
        legend.position = "bottom"
      )
ggsave(paste(global_var$global$path_microglai_statistics, "fraction_replicates_seperated.png", sep = "/"), p , width = 3.5 , height = 5 , units = "in")

############ Box Plot for all microglia
############ generate meta data, for statistical testing and box plot 
integrated.meta.stat <- integrated.strain@meta.data %>%
                mutate( strain=factor(global_var$global$strain, levels = c("WT", "AZT")),
                        new_clusters=factor(seurat_clusters, levels = c("1","2","3","4","5","6", "7"))) %>%
                group_by(strain, new_clusters) %>%
                summarise(Med_nFeature=median(nFeature_RNA),
                          Med_percent_mt= median(percent.mt),
                          med_percent.microglia=median(percent.microglia),
                          
                          N=n()) %>%
  group_by(Percent= N/sum(N)*100)

integratedPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#f0E442", "#0072B2", "#D55E00")  

p <- integrated.meta.stat %>%
  ggplot(aes(y=Percent, x= strain, color = strain )) +
  geom_boxplot(outlier.size = 0, alph= 0.5) +
  geom_point(aes(color = strain), position = position_jitterdodge(), alpha=0.8) +
  scale_color_manual(values = integratedPalette) +
  theme_bw() +
  facet_grid(new_clusters ~ strain, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size= 10),
        strip.text = element_text( face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom"
        )

ggsave(paste(global_var$global$path_microglai_statistics, "cluster_box_all.png", sep = "/"), p , width = 4, height = 10, units = "in")


####### Perform two-way ANOVA to determine the effect of strain on the percent of microglia subclusters

clusters <- unique(integrated.meta.stat$new_clusters) %>% as.list()
data  = integrated.meta.stat %>%  filter(new_clusters %in% clusters[[1]])
aov_object <-aov( integrated.meta.stat$Percent ~ integrated.meta.stat$new_clusters*integrated.meta.stat$med_percent.microglia, data = data)
aov.pvals <- summary(aov_object)
aov.pvals <- aov.pvals %>% t()


######### check the statistics on nFeature, percent of microglia, and percent of ribosomal genes for each cluster 

#### nFeature
aov_stat <- aov(Med_nFeature ~ new_clusters, data = integrated.meta.stat)
aov_table <-TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>%
  mutate(comparison = paste(" ", comparison , sep = ""), Significance=ifelse(p.adj <0.05, "S", "NS"))
write_delim(aov_table, path = paste(global_var$global$path_microglai_statistics, "Med_nFeature_comp_cluster.txt", sep = "/"), delim = "\t")
filter(aov_table, Significance == "S")


###### percent of mitochodria
aov_stat = aov(Med_percent_mt ~ new_clusters, data = integrated.meta.stat)
aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data_frame() %>% rownames_to_column(var = "comparison") %>%
  mutate(comparison = paste(" ", comparison, sep = " "), Significance=ifelse("p.adj" < 0.05 , "S", "NS"))
filter(aov_table, Significance =="S")

#write_delim(aov_table, path = paste(global_var$global$path_microglai_statistics, "Med_percent_mt_comp_cluster.txt", sep = "/"), delim = "\t")

################## percent of Microglia
aov_stat = aov(med_percent.microglia ~ new_clusters, data = integrated.meta.stat)

aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var =  "comparison ") %>%
  mutate(comparison = paste(" ", comparison, sep = " "), Significance=ifelse("p.adj" < 0.05 , "S", "NS"))
