library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(tidyverse)
library(patchwork)
library(viridis)
library(scCustomize)

DefaultAssay(all)="RNA"

all<- all%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




all<- all%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

all<- all%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()





FeaturePlot(all, features =c("EPCAM", "PTPRC", "PECAM1", "CLDN5", "COL1A1", "ACTA2"), cols=c("lightgrey", "red"), min.cutoff = "q10")



endo=WhichCells(all, idents=c(8, 16))


fib=WhichCells(all, idents=11)


smc=WhichCells(all, idents=15)


DimPlot(all, label=T)


all=SetIdent(all, cells=endo, value="Endothelial")
all=SetIdent(all, cells=fib, value="Fibroblast")
all=SetIdent(all, cells=smc, value="SMC")


DimPlot(all, label=T)



saveRDS(all, file="allsymbolupdate2.rds" )


epim=subset(all, idents=c("Endothelial", "Fibroblast", "SMC"), invert=T)


epim<- epim%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




epim<- epim%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)



epim<- epim%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

DimPlot(epim, label=T)



saveRDS(epim, file="epimupdate.rds")
FeaturePlot(LongCovid, features=c("LYZ","FCGR3A","CD14", "FCER1A","LILRA4","FCGR3B", "CD3E", "CD3D", "CD4", "CD8A" ,"CD8B", "GNLY","TYROBP","NCAM1","KLRD1", "CD19", "MS4A1", "PPBP"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=4)


FeaturePlot(LongCovid, features=c( "GNLY", "NKG7","GZMA", "GZMB", "GZMH"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=3)

FeaturePlot(LongCovid, features =c("CPA3", "KIT"), cols=c("lightgrey", "red"), min.cutoff = "q10")




FeaturePlot(LongCovid, features =c("FCER1A", "CD1C"), cols=c("lightgrey", "red"), min.cutoff = "q10")


FeaturePlot(LongCovid, features =c("LILRA4", "IL3RA"), cols=c("lightgrey", "red"), min.cutoff = "q10")



FeaturePlot(epim, features =c( "PTPRC", "EPCAM"), cols=c("lightgrey", "red"), min.cutoff = "q10")


FeaturePlot(epim, features =c("SFTPC", "SFTPB","SCGB1A1", "SCGB3A1","SCGB3A2", "MUC5B", "MUC5AC","LTF","PRB3", "FOXJ1","CAPS", "C9orf24", "MCIDAS", "KRT5", "KRT15", "KRT17", "AGER", "CAV1", "FOXI1"), cols=c("lightgrey", "red"), min.cutoff = "q10")




ep=subset(epim, idents=c(0, 16, 22, 25, 19, 24, 3, 10, 7, 8, 1, 18, 13, 26, 21, 12))
DimPlot(ep)
rm(all, epim, epim)



ep<- ep%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




ep<- ep%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)



ep<- ep%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()


DimPlot(ep, label=T)


FeaturePlot(ep, features =c("SFTPC", "SFTPB","SCGB1A1", "SCGB3A2", "MUC5B", "MUC5AC","LTF","PRB3", "FOXJ1","CAPS", "C9orf24", "KRT5", "KRT15", "KRT17", "AGER", "CAV1"), cols=c("lightgrey", "red"), order=T, min.cutoff = "q10")



FeaturePlot(ep, features=c("LTF", "PRB3"), min.cutoff = "q10", max.cutoff = "q95", order=T)

FeaturePlot(ep, features=c("EPCAM", "TOP2A", "PTPRC", "MKI67"), min.cutoff = "q10", max.cutoff = "q95", order=T)

table(ep$celltype)

c21=WhichCells(ep, idents=21)
FeaturePlot(ep, features=c("CD19", "MS4A1", "JCHAIN", "LYZ"), min.cutoff = "q10", max.cutoff = "q95", order=T)



ep=RenameIdents(ep, "0"="AT2", "4"="AT2", "14"="AT2", "16"="AT2", "2"="SCGB3A2+", "8"="Club", "6"="Goblet", "12"="KRT17+ KRT15- Basal", "7"="KRT17+ KRT15+ Basal",
                "15"="Transitional", "9"="AT1", "18"="pre-Ciliated", "21"="Serous", "1"="Ciliated", "5"="Ciliated", "3"="Ciliated", "11"="Ciliated", "17"="Ciliated", "19"="Ciliated", "20"="Ciliated", "22"="Cycling", "13"="Cycling", "10"="Cycling"   )


levels(ep)=c("AT1", "AT2", "Ciliated", "Club", "Cycling", "Goblet", "KRT17+ KRT15- Basal", "KRT17+ KRT15+ Basal", "pre-Ciliated", "SCGB3A2+", "Serous", "Transitional")



ep$celltype=Idents(ep)

saveRDS(ep, file="epudate.rds")

table(ep$celltype2)

ep<- ep%>% 
  FindClusters(resolution = 1.0) %>% 
  identity()




cyc=subset(ep, idents="Cycling")





cyc<- cyc%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

DimPlot(cyc, label=T)

FeaturePlot(cyc, features =c("SFTPC", "SFTPB","SCGB1A1", "SCGB3A2", "MUC5B", "MUC5AC","LTF","PRB3", "FOXJ1","CAPS", "C9orf24", "KRT5", "KRT15", "KRT17", "AGER", "CAV1"), cols=c("lightgrey", "red"), order=T, min.cutoff = "q10")

FeaturePlot(cyc, features=c("LYZ","FCGR3A","CD14", "FCER1A","LILRA4","FCGR3B", "CD3E", "CD3D", "CD4", "CD8A" ,"CD8B", "GNLY","TYROBP","NCAM1","KLRD1", "CD19", "MS4A1", "PPBP", "CPA3", "KIT"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=4)


FeaturePlot(cyc, features=c("LTF", "PRB3"), min.cutoff = "q10", max.cutoff = "q95", order=T)

FeaturePlot(cyc, features=c("EPCAM", "TOP2A", "PTPRC", "MKI67", "TUBB", "TUBA1B"), min.cutoff = "q10", max.cutoff = "q95", order=T,cols=c("lightgrey", "red"))

cyc<- cyc%>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()



DimPlot(cyc, label=T)

DimPlot(ep, split.by="diagnosis")

im_ep=WhichCells(cyc, idents=c(1,  7))


at2=WhichCells(cyc, idents=c(4))
mac=WhichCells(cyc, idents=c(6))

DimPlot(ep, cells.highlight = im_ep)


DimPlot(ep, cells.highlight = at2)

ep=SetIdent(ep, cells=at2, value="AT2")

ep$celltype2=Idents(ep)

FeaturePlot(ep, features=c("EPCAM", "TOP2A", "PTPRC", "MKI67"), min.cutoff = "q10", max.cutoff = "q95", order=T,cols=c("lightgrey", "red"))


DimPlot()

DimPlot(epim, cells.highlight = mac)


ep=SetIdent(ep, cells=mac, value="Macrophage")

ep$celltype2=Idents(ep)

DimPlot(ep)

ep=SetIdent(ep, cells=im_ep, value="im-ep")

ep$celltype2=Idents(ep)


cells=WhichCells(ep)
im=subset(epim, cells=cells, invert=T)

















DimPlot(im, label=T)


im<- im%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




im<- im%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)



im<- im%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

DimPlot(im, label=T)






FeaturePlot(im, features=c("EPCAM", "TOP2A", "PTPRC", "MKI67"), min.cutoff = "q10", max.cutoff = "q95", order=T)



FeaturePlot(im, features=c("CD1D", "FCER1A"), min.cutoff = "q10", max.cutoff = "q95", order=T)





FeaturePlot(im, features=c("LYZ","FCGR3A","CD14", "FCER1A","LILRA4","FCGR3B", "CD3E", "CD3D", "CD4", "CD8A" ,"CD8B", "GNLY","TYROBP","NCAM1","KLRD1", "CD19", "MS4A1", "PPBP", "CPA3", "KIT"), cols=c("lightgrey", "red"), max.cutoff = "q90", ncol=4)




FeaturePlot(im, features =c("SFTPC", "SFTPB","SCGB1A1", "SCGB3A2", "MUC5B", "MUC5AC","LTF","PRB3", "FOXJ1","CAPS", "C9orf24", "KRT5", "KRT15", "KRT17", "AGER", "CAV1"), cols=c("lightgrey", "red"),  min.cutoff = "q10")






FeaturePlot(im, features =c("EPCAM", "PTPRC"), cols=c("lightgrey", "red"),  min.cutoff = "q10")





im=RenameIdents(im, "0"="Macrophage", "1"="Macrophage", "2"="Macrophage", "19"="Macrophage", "16"="Macrophage", "3"="Macrophage", "9"="Macrophage", "13"="Macrophage", "6"="cDC", "17"="cDC", "8"="Neutrophil", "15"="Plasma", "12"="B cell", 
                "7"="NK cell", "4"="CD8 T cell", "5"="CD4 T cell", "14"="CD4 T cell", "10"="Mast cell", "11"="im-ep", "18"="pDC")
DimPlot(im)




levels(im)=c("B cell", "CD4 T cell", "CD8 T cell", "cDC",  "im-ep" , "Mast cell", "Macrophage",   "Neutrophil", "NK cell", "pDC",  "Plasma" )




im$celltype2=Idents(im)


saveRDS(im, file="imupdate.rds")





epmeta=ep@meta.data
immeta=im@meta.data

celltype2=epmeta[ , "celltype2", drop=F]




imcelltype2=immeta[ , "celltype2", drop=F]









cell=rbind(celltype2, imcelltype2)


epim=AddMetaData(epim, cell, col.name="celltype2")



DimPlot(epim, group.by="celltype2", label=T)

DimPlot(cyc, label=T)


Idents(epim)="celltype2"

im_ep=WhichCells(epim, idents="im-ep")
DimPlot(epim, cells.highlight = cycling)


cycling=WhichCells(ep, idents="Cycling")



ep=subset(ep, idents="Macrophage", invert=T)




levels(ep)=c("AT1", "AT2", "Ciliated", "Club", "Cycling", "Goblet", "im-ep","KRT17+ KRT15- Basal", "KRT17+ KRT15+ Basal", "pre-Ciliated", "SCGB3A2+", "Serous", "Transitional")
ep$celltype2=Idents(ep)
saveRDS(ep, file="epupdate.rds")


saveRDS(im, file="imupdate2.rds")



imm=subset(epim, idents=c("B cell", "CD4 T cell", "CD8 T cell", "cDC", "Cycling Immune", "im-ep" , "Mast cell", "Macrophage",   "Neutrophil", "NK cell", "pDC",  "Plasma"))



imm<- imm%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




imm<- imm%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)



imm<- imm%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()


saveRDS(imm, file="imm.rds")





DimPlot(imm, label=T, group.by="celltype3")
DimPlot(imm, label=T)
imm$RNA_snn_res.0.6
Idents(imm)="RNA_snn_res.0.6"

c11=WhichCells(imm, idents=11)
Idents(imm)=imm$celltype3

mac=subset(imm, idents="Macrophage")




DimPlot(ep)
ep$celltype2=Idents(ep)

epmeta=epudate@meta.data
immeta=im@meta.data
table(ep$celltype2)
table(im$celltype2)

table(epim$celltype2)

epcelltype2=epmeta[, "celltype2", drop=F]
imcelltype2=immeta[, "celltype2", drop=F]

celltype2=rbind(epcelltype2, imcelltype2)


epimep=WhichCells(ep, idents="im-ep")
epudate=SetIdent(epudate, cells=epimep, value="im-ep")
epudate=SetIdent(epudate, cells=mac, value="Macrophage")

epudate$celltype3=Idents(epudate)
im$celltype3=Idents(im)

epmeta=epudate@meta.data
immeta=im@meta.data




epcelltype3=epmeta[, "celltype3", drop=F]
imcelltype3=immeta[, "celltype3", drop=F]

celltype3=rbind(epcelltype3, imcelltype3)


epim=AddMetaData(epim, celltype3, col.name = "celltype3")

table(epim$celltype3)



DimPlot(epim, label=T, group.by="celltype3")
imep=WhichCells(epim, idents="im-ep")


DimPlot(epim, cells.highlight = imep)

Idents(epim)=epim$celltype3

cycling=subset(epim, idents="Cycling")




levels(epudate)=c("AT1", "AT2", "Ciliated", "Club", "Cycling Epithelial", "Cycling Immune", "Goblet", "im-ep","KRT17+ KRT15- Basal", "KRT17+ KRT15+ Basal", "pre-Ciliated", "SCGB3A2+", "Serous", "Transitional", "Macrophage")
epudate$celltype3=Idents(epudate)

saveRDS(epudate, file="epupdate.rds")


Idents(epim)="celltype3"

levels(epim)=c("AT1", "AT2", "Ciliated", "Club", "Cycling Epithelial",  "Goblet", "KRT17+ KRT15- Basal", "KRT17+ KRT15+ Basal", "pre-Ciliated", "SCGB3A2+", "Serous", "Transitional", "B cell", "CD4 T cell", "CD8 T cell", "cDC", "Cycling Immune", "im-ep" , "Mast cell", "Macrophage",   "Neutrophil", "NK cell", "pDC",  "Plasma")

epim$celltype3=Idents(epim)
saveRDS(epim, file="epim.rds")


ep=subset(epudate, idents=c("Cycling Immune",  "im-ep", "Macrophage" ), invert=T)



ep<- ep%>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = all@var.genes, npcs = 30, verbose = T)
###########




ep<- ep%>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)



ep<- ep%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

DimPlot(ep, label=T, group.by="celltype3")
DimPlot(ep, label=T)

saveRDS(ep, file="epupdate2.rds")


epcells=WhichCells(ep)


imm=subset(epim, cells=ep, invert=T)


imm<- imm%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()
DimPlot(imm, label=T)
Idents(imm)=imm$celltype3



DimPlot(imm)
saveRDS(imm, file="imm2.rds")

epim=RenameIdents(epim, "KRT17+ KRT15+ Basal"="Basal",  "KRT17+ KRT15- Basal"="Basal", "pre-Ciliated"="Ciliated")


levels(epim)=c("AT1",    "AT2", "Basal","Ciliated", "Club",   "Cycling Epithelial", "Goblet",  "SCGB3A2+" , "Serous", "Transitional", 
               "B cell" ,  "CD4 T cell",   "CD8 T cell" , "cDC",  "Cycling Immune", "ImEp", "Macrophage", "Mast cell" ,          "Neutrophil",          "NK cell", "pDC" ,                "Plasma"  )
epim$celltype5=Idents(epim)


saveRDS(epim, file="epim2.rds")

library(future)
options(future.globals.maxSize= 2000000000)
plan("multiprocess", workers = 4)

epmarkers=FindAllMarkers(ep, only.pos = T)
saveRDS(epmarkers, file="epmarkers.rds")

top10 <- epmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- epmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DoHeatmap(subset(ep, downsample = 1000), features = top20$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))

DoHeatmap(subset(ep, downsample = 1000), features = top10$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))


top10_markers <- Extract_Top_Markers(marker_dataframe = epmarkers, num_genes = 10, named_vector = FALSE,
                                    make_unique = TRUE)

Clustered_DotPlot(seurat_object = ep, features = top10_markers)



library(future)
options(future.globals.maxSize= 2000000000)
plan("multiprocess", workers = 4)

immmarkers=FindAllMarkers(imm, only.pos = T)
saveRDS(immmarkers, file="immmarkers.rds")

top10 <- immmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- immmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DoHeatmap(subset(imm, downsample = 500), features = top20$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))

DoHeatmap(subset(imm, downsample = 500), features = top10$gene) +  scale_fill_gradientn(colors = c("lightblue", "white", "red"))


top10_markers <- Extract_Top_Markers(marker_dataframe = immmarkers, num_genes = 10, named_vector = FALSE,
                                     make_unique = TRUE)

Clustered_DotPlot(seurat_object = imm, features = top10_markers)




















