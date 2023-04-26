
library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
DimPlot(imm)


Mac=subset(imm, idents="Macrophage")


Mac<- Mac%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = Mac@var.genes, npcs = 30, verbose = T)
###########





Mac<- Mac%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

Mac<- Mac%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


DimPlot(Mac, label=T)

FeaturePlot(Mac, features=c("LYZ","FCGR3A","CD14", "FCER1A","CD1C","FCGR3B", "CD3E", "CXCR1", "CX3CR1", "CCR2", "C1QC", "CD68"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=4)
FeaturePlot(Mac, features=c("MAFB", "MAF", "ITGAM", "MERTK"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=4)


FeaturePlot(Mac, features=c( "GNLY", "NKG7","GZMA", "GZMB", "GZMH"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=3)


FeaturePlot(Mac, features=c( "FCER1A", "CD1C"), cols=c("lightgrey", "red"), max.cutoff = "q90")


FeaturePlot(Mac, features=c( "FABP4", "C1QB", "ITGAM", "MERTK", "SPP1", "FCN1"), cols=c("lightgrey", "red"), max.cutoff = "q90")


FeaturePlot(Mac, features=c( "FABP4", "APOC1", "MARCO", "MERTK", "SPP1", "TREM2", "FCN1", "S100A8", "CCL2", "CCL3", "FCGR3A", "CD14"), cols=c("lightgrey", "red"), max.cutoff = "q90")

FeaturePlot(Mac, features=c( "FABP4", "APOC1", "MARCO", "MERTK", "SPP1", "TREM2", "FCN1", "S100A8", "CCL2", "CCL3", "FCGR3A", "CD14"), cols=c("lightgrey", "red"), max.cutoff = "q95", min.cutoff = "q10", order=T)






FeaturePlot(Mac, features=c( "CD163"), cols=c("lightgrey", "red"), max.cutoff = "q90")





Mac=RenameIdents(Mac, "0"="FABP4+ Mac", "4"="FABP4+ Mac", "6"="FABP4+ Mac", "10"="FABP4+ Mac", "1"="FCN1+ CD14+ Mac", "8"="FCN1+ CD16+ Mac", '3'="Intermediate Mac",
                '9'="Intermediate Mac", '7'="Intermediate Mac", "5"="MERTK+ Mac", "2"="SPP1+ Mac" , "11"="T monocyte complex")
Mac$celltype=Idents(Mac)


Idents(Mac)=Mac$RNA_snn_res.0.4

Mac=RenameIdents(Mac, "0"="FABP4+ Mac", "4"="FABP4+ Mac", "6"="FABP4+ Mac", "10"="FABP4+ Mac", "1"="FCN1+ CD14+ Mac", "8"="FCN1+ CD16+ Mac", '3'="Intermediate Mac",
                '9'="Intermediate Mac", '7'="Intermediate Mac", "5"="MERTK+ SPP1+ Mac", "2"="MERTK+ SPP1+ Mac" , "11"="T monocyte complex")


Mac$celltype2=Idents(Mac)
saveRDS(Mac, file="Mac.rds")

DimPlot(Mac)









library(future)
options(future.globals.maxSize= 2000000000)
plan("multiprocess", workers = 4)

Macmarkers=FindAllMarkers(Mac, only.pos = T)

top10 <- Macmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(subset(Mac, downsample = 1000), features = top10$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))

top20 <- Macmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DoHeatmap(subset(Mac, downsample = 1000), features = top20$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))


saveRDS(Macmarkers, file="Macmarkers.rds")

write.csv(Macmarkers2, "Macmarkers.csv")

DDimPlot(Mac, split.by="diagnosis")



saveRDS(Mac, file="Mac.rds")

Mac=subset(Mac, idents="T monocyte complex", invert=T)

levels(Mac)
levels(Mac)=c("FABP4+ Mac" ,   "FCN1+ CD16+ Mac", "FCN1+ CD14+ Mac", "Intermediate Mac", "MERTK+ Mac", "SPP1+ Mac"  )
library(future)
options(future.globals.maxSize= 2000000000)
plan("multiprocess", workers = 4)
Idents(Mac)=Mac$celltype4
Macmarkers2=FindAllMarkers(Mac, only.pos = T)

top10 <- Macmarkers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(subset(Mac, downsample = 1000), features = top10$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))
saveRDS(Macmarkers, file="Macmarkers2.rds")


saveRDS(Mac, file="Mac2.rds")





library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)

top10_markers <- Extract_Top_Markers(marker_dataframe = Macmarkers, num_genes = 10, named_vector = FALSE,
                                     make_unique = TRUE)

Clustered_DotPlot(seurat_object = Mac, features = top10_markers)




