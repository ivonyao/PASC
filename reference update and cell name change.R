library(stringr)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(AnnotationHub)
library(ensembldb)
library(cowplot)
library(ggplot2)
library(scales)
library(Seurat.utils)

ldat7 <- ReadVelocity(file = "h7_2.loom")

colnames(ldat7[["spliced"]])=str_replace(colnames(ldat7[["spliced"]]), "h7_2:", "lc1a_")
colnames(ldat7[["spliced"]])=str_replace(colnames(ldat7[["spliced"]]), "x", "-1")
colnames(ldat7[["spliced"]])



colnames(ldat7[["unspliced"]])=str_replace(colnames(ldat7[["unspliced"]]), "h7_2:", "lc1a_")
colnames(ldat7[["unspliced"]])=str_replace(colnames(ldat7[["unspliced"]]), "x", "-1")
colnames(ldat7[["unspliced"]])

colnames(ldat7[["ambiguous"]])=str_replace(colnames(ldat7[["ambiguous"]]), "h7_2:", "lc1a_")
colnames(ldat7[["ambiguous"]])=str_replace(colnames(ldat7[["ambiguous"]]), "x", "-1")
colnames(ldat7[["ambiguous"]])

h7=as.Seurat(x = ldat7)


h7_2=Read10X_h5('h7_2_filtered_feature_bc_matrix.h5')



colnames(h7_2)=paste("lc1a", colnames(h7_2), sep="_")



h7_2=CreateSeuratObject(counts=h7_2, project="h7_2")

h7[["RNA"]] <- h7_2[["RNA"]]
DefaultAssay(h7)="RNA"










h7[["percent.mt"]] <- PercentageFeatureSet(h7, pattern = "^MT-")
VlnPlot(h7,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h7))
q = (100 - perc)/100
feature_limit <- quantile(h7$nFeature_RNA, q)

h7 <- subset(h7, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

h7$patient="LC1"


h7$diagnosis="LongCovid"

h7$region="apical"
saveRDS(h7, file="LC1apical.rds")





ldat8<- ReadVelocity(file = "h8_2.loom")


colnames(ldat8[["spliced"]])=str_replace(colnames(ldat8[["spliced"]]), "h8_2:", "lc1b_")
colnames(ldat8[["spliced"]])=str_replace(colnames(ldat8[["spliced"]]), "x", "-1")
colnames(ldat8[["spliced"]])



colnames(ldat8[["unspliced"]])=str_replace(colnames(ldat8[["unspliced"]]), "h8_2:", "lc1b_")
colnames(ldat8[["unspliced"]])=str_replace(colnames(ldat8[["unspliced"]]), "x", "-1")
colnames(ldat8[["unspliced"]])

colnames(ldat8[["ambiguous"]])=str_replace(colnames(ldat8[["ambiguous"]]), "h8_2:", "lc1b_")
colnames(ldat8[["ambiguous"]])=str_replace(colnames(ldat8[["ambiguous"]]), "x", "-1")
colnames(ldat8[["ambiguous"]])


h8=as.Seurat(x = ldat8)

h8_2=Read10X_h5('h8_2_filtered_feature_bc_matrix.h5')



colnames(h8_2)=paste("lc1b", colnames(h8_2), sep="_")



h8_2=CreateSeuratObject(counts=h8_2, project="h8_2")


h8[["RNA"]] <- h8_2[["RNA"]]
DefaultAssay(h8)="RNA"



h8[["percent.mt"]] <- PercentageFeatureSet(h8, pattern = "^MT-")
VlnPlot(h8,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h8))
q = (100 - perc)/100
feature_limit <- quantile(h8$nFeature_RNA, q)

h8 <- subset(h8, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

h8$patient="LC1"


h8$diagnosis="LongCovid"

h8$region="basal"
saveRDS(h8, file="LC1basal.rds")


ldat9<- ReadVelocity(file = "h9_2.loom")

colnames(ldat9[["spliced"]])=str_replace(colnames(ldat9[["spliced"]]), "h9_2:", "lc2a_")
colnames(ldat9[["spliced"]])=str_replace(colnames(ldat9[["spliced"]]), "x", "-1")
colnames(ldat9[["spliced"]])



colnames(ldat9[["unspliced"]])=str_replace(colnames(ldat9[["unspliced"]]), "h9_2:", "lc2a_")
colnames(ldat9[["unspliced"]])=str_replace(colnames(ldat9[["unspliced"]]), "x", "-1")
colnames(ldat9[["unspliced"]])

colnames(ldat9[["ambiguous"]])=str_replace(colnames(ldat9[["ambiguous"]]), "h9_2:", "lc2a_")
colnames(ldat9[["ambiguous"]])=str_replace(colnames(ldat9[["ambiguous"]]), "x", "-1")
colnames(ldat9[["ambiguous"]])

h9_2=Read10X_h5('h9_2_filtered_feature_bc_matrix.h5')



colnames(h9_2)=paste("lc2a", colnames(h9_2), sep="_")



h9_2=CreateSeuratObject(counts=h9_2, project="h9_2")

h9=as.Seurat(x = ldat9)

h9[["RNA"]] <- h9_2[["RNA"]]
DefaultAssay(h9)="RNA"

h9[["percent.mt"]] <- PercentageFeatureSet(h9, pattern = "^MT-")
VlnPlot(h9,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h9))
q = (100 - perc)/100
feature_limit <- quantile(h9$nFeature_RNA, q)

h9 <- subset(h9, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

h9$patient="LC2"


h9$diagnosis="LongCovid"

h9$region="apical"
saveRDS(h9, file="LC2apical.rds")



ldat10 <- ReadVelocity(file = "h10_2.loom")


colnames(ldat10[["spliced"]])=str_replace(colnames(ldat10[["spliced"]]), "h10_2:", "lc2b_")
colnames(ldat10[["spliced"]])=str_replace(colnames(ldat10[["spliced"]]), "x", "-1")
colnames(ldat10[["spliced"]])



colnames(ldat10[["unspliced"]])=str_replace(colnames(ldat10[["unspliced"]]), "h10_2:", "lc2b_")
colnames(ldat10[["unspliced"]])=str_replace(colnames(ldat10[["unspliced"]]), "x", "-1")
colnames(ldat10[["unspliced"]])

colnames(ldat10[["ambiguous"]])=str_replace(colnames(ldat10[["ambiguous"]]), "h10_2:", "lc2b_")
colnames(ldat10[["ambiguous"]])=str_replace(colnames(ldat10[["ambiguous"]]), "x", "-1")
colnames(ldat10[["ambiguous"]])

h10=as.Seurat(x = ldat10)
h10_2=Read10X_h5('h10_2_filtered_feature_bc_matrix.h5')



colnames(h10_2)=paste("lc2b", colnames(h10_2), sep="_")



h10_2=CreateSeuratObject(counts=h10_2, project="h10_2")

h10[["RNA"]] <- h10_2[["RNA"]]
DefaultAssay(h10)="RNA"


h10[["percent.mt"]] <- PercentageFeatureSet(h10, pattern = "^MT-")

VlnPlot(h10,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h10))
q = (100 - perc)/100
feature_limit <- quantile(h10$nFeature_RNA, q)

h10 <- subset(h10, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

h10$patient="LC2"


h10$diagnosis="LongCovid"

h10$region="basal"
saveRDS(h10, file="LC2basal_velo.rds")




ldat11<- ReadVelocity(file = "h11_2.loom")

ldat11 <- ReadVelocity(file = "h11_2.loom")

colnames(ldat11[["spliced"]])=str_replace(colnames(ldat11[["spliced"]]), "h11_2:", "lc3a_")
colnames(ldat11[["spliced"]])=str_replace(colnames(ldat11[["spliced"]]), "x", "-1")
colnames(ldat11[["spliced"]])



colnames(ldat11[["unspliced"]])=str_replace(colnames(ldat11[["unspliced"]]), "h11_2:", "lc3a_")
colnames(ldat11[["unspliced"]])=str_replace(colnames(ldat11[["unspliced"]]), "x", "-1")
colnames(ldat11[["unspliced"]])

colnames(ldat11[["ambiguous"]])=str_replace(colnames(ldat11[["ambiguous"]]), "h11_2:", "lc3a_")
colnames(ldat11[["ambiguous"]])=str_replace(colnames(ldat11[["ambiguous"]]), "x", "-1")
colnames(ldat11[["ambiguous"]])

h11=as.Seurat(x = ldat11)
h11_2=Read10X_h5('h11_2_filtered_feature_bc_matrix.h5')



colnames(h11_2)=paste("lc3a", colnames(h11_2), sep="_")



h11_2=CreateSeuratObject(counts=h11_2, project="h11_2")

h11[["RNA"]] <- h11_2[["RNA"]]
DefaultAssay(h11)="RNA"



h11[["percent.mt"]] <- PercentageFeatureSet(h11, pattern = "^MT-")
VlnPlot(h11,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h11))
q = (100 - perc)/100
feature_limit <- quantile(h11$nFeature_RNA, q)

h11 <- subset(h11, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

h11$patient="LC3"


h11$diagnosis="LongCovid"

h11$region="apical"
saveRDS(h11, file="LC3apical.rds")



ldat12<- ReadVelocity(file = "h12_2.loom")


ldat12 <- ReadVelocity(file = "h12_2.loom")

colnames(ldat12[["spliced"]])=str_replace(colnames(ldat12[["spliced"]]), "h12_2:", "lc3b_")
colnames(ldat12[["spliced"]])=str_replace(colnames(ldat12[["spliced"]]), "x", "-1")
colnames(ldat12[["spliced"]])



colnames(ldat12[["unspliced"]])=str_replace(colnames(ldat12[["unspliced"]]), "h12_2:", "lc3b_")
colnames(ldat12[["unspliced"]])=str_replace(colnames(ldat12[["unspliced"]]), "x", "-1")
colnames(ldat12[["unspliced"]])

colnames(ldat12[["ambiguous"]])=str_replace(colnames(ldat12[["ambiguous"]]), "h12_2:", "lc3b_")
colnames(ldat12[["ambiguous"]])=str_replace(colnames(ldat12[["ambiguous"]]), "x", "-1")
colnames(ldat12[["ambiguous"]])

h12=as.Seurat(x = ldat12)
h12_2=Read10X_h5('h12_2_filtered_feature_bc_matrix.h5')



colnames(h12_2)=paste("lc3b", colnames(h12_2), sep="_")



h12_2=CreateSeuratObject(counts=h12_2, project="h12_2")

h12[["RNA"]] <- h12_2[["RNA"]]
DefaultAssay(h12)="RNA"




h12[["percent.mt"]] <- PercentageFeatureSet(h12, pattern = "^MT-")
VlnPlot(h12,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(h12))
q = (100 - perc)/100
feature_limit <- quantile(h12$nFeature_RNA, q)

h12 <- subset(h12, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)



h12$patient="LC3"


h12$diagnosis="LongCovid"

h12$region="basal"

saveRDS(h12, file="LC3basal.rds")



ldatlc4<- ReadVelocity(file = "lc4_2.loom")


ldatlc4 <- ReadVelocity(file = "lc4_2.loom")

colnames(ldatlc4[["spliced"]])=str_replace(colnames(ldatlc4[["spliced"]]), "lc4_2:", "lc4_")
colnames(ldatlc4[["spliced"]])=str_replace(colnames(ldatlc4[["spliced"]]), "x", "-1")
colnames(ldatlc4[["spliced"]])



colnames(ldatlc4[["unspliced"]])=str_replace(colnames(ldatlc4[["unspliced"]]), "lc4_2:", "lc4_")
colnames(ldatlc4[["unspliced"]])=str_replace(colnames(ldatlc4[["unspliced"]]), "x", "-1")
colnames(ldatlc4[["unspliced"]])

colnames(ldatlc4[["ambiguous"]])=str_replace(colnames(ldatlc4[["ambiguous"]]), "lc4_2:", "lc4_")
colnames(ldatlc4[["ambiguous"]])=str_replace(colnames(ldatlc4[["ambiguous"]]), "x", "-1")
colnames(ldatlc4[["ambiguous"]])

lc4=as.Seurat(x = ldatlc4)
lc4_2=Read10X_h5('lc4_2_filtered_feature_bc_matrix.h5')



colnames(lc4_2)=paste("lc4", colnames(lc4_2), sep="_")



lc4_2=CreateSeuratObject(counts=lc4_2, project="lc4_2")

lc4[["RNA"]] <- lc4_2[["RNA"]]
DefaultAssay(lc4)="RNA"


lc4[["percent.mt"]] <- PercentageFeatureSet(lc4, pattern = "^MT-")
VlnPlot(lc4,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(lc4))
q = (100 - perc)/100
feature_limit <- quantile(lc4$nFeature_RNA, q)

lc4 <- subset(lc4, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)


lc4$patient="LC4"


lc4$diagnosis="LongCovid"

lc4$region="apical_basal"
saveRDS(lc4, file="LC4.rds")


ldatlc5a <- ReadVelocity(file = "lc5apical.loom")
colnames(ldatlc5a[["spliced"]])
colnames(ldatlc5a[["spliced"]])=str_replace(colnames(ldatlc5a[["spliced"]]), "lc5apical:", "lc5a_")
colnames(ldatlc5a[["spliced"]])=str_replace(colnames(ldatlc5a[["spliced"]]), "x", "-1")
colnames(ldatlc5a[["spliced"]])



colnames(ldatlc5a[["unspliced"]])=str_replace(colnames(ldatlc5a[["unspliced"]]), "lc5apical:", "lc5a_")
colnames(ldatlc5a[["unspliced"]])=str_replace(colnames(ldatlc5a[["unspliced"]]), "x", "-1")
colnames(ldatlc5a[["unspliced"]])

colnames(ldatlc5a[["ambiguous"]])=str_replace(colnames(ldatlc5a[["ambiguous"]]), "lc5apical:", "lc5a_")
colnames(ldatlc5a[["ambiguous"]])=str_replace(colnames(ldatlc5a[["ambiguous"]]), "x", "-1")
colnames(ldatlc5a[["ambiguous"]])

lc5a=as.Seurat(x = ldatlc5a)


lc5a_2=Read10X_h5('lc5a_filtered_feature_bc_matrix.h5')



colnames(lc5a_2)=paste("lc5a", colnames(lc5a_2), sep="_")



lc5a_2=CreateSeuratObject(counts=lc5a_2, project="lc5a_2")

lc5a[["RNA"]] <- lc5a_2[["RNA"]]
DefaultAssay(lc5a)="RNA"










lc5a[["percent.mt"]] <- PercentageFeatureSet(lc5a, pattern = "^MT-")
VlnPlot(lc5a,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(lc5a))
q = (100 - perc)/100
feature_limit <- quantile(lc5a$nFeature_RNA, q)

lc5a <- subset(lc5a, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

lc5a$patient="LC5"


lc5a$diagnosis="LongCovid"

lc5a$region="apical"
saveRDS(lc5a, file="LC5apical.rds")


ldatlc5b <- ReadVelocity(file = "lc5basal.loom")
colnames(ldatlc5b[["spliced"]])
colnames(ldatlc5b[["spliced"]])=str_replace(colnames(ldatlc5b[["spliced"]]), "lc5basal:", "lc5b_")
colnames(ldatlc5b[["spliced"]])=str_replace(colnames(ldatlc5b[["spliced"]]), "x", "-1")
colnames(ldatlc5b[["spliced"]])



colnames(ldatlc5b[["unspliced"]])=str_replace(colnames(ldatlc5b[["unspliced"]]), "lc5basal:", "lc5b_")
colnames(ldatlc5b[["unspliced"]])=str_replace(colnames(ldatlc5b[["unspliced"]]), "x", "-1")
colnames(ldatlc5b[["unspliced"]])

colnames(ldatlc5b[["ambiguous"]])=str_replace(colnames(ldatlc5b[["ambiguous"]]), "lc5basal:", "lc5b_")
colnames(ldatlc5b[["ambiguous"]])=str_replace(colnames(ldatlc5b[["ambiguous"]]), "x", "-1")
colnames(ldatlc5b[["ambiguous"]])

lc5b=as.Seurat(x = ldatlc5b)


lc5b_2=Read10X_h5('lc5b_filtered_feature_bc_matrix.h5')



colnames(lc5b_2)=paste("lc5b", colnames(lc5b_2), sep="_")



lc5b_2=CreateSeuratObject(counts=lc5b_2, project="lc5b_2")

lc5b[["RNA"]] <- lc5b_2[["RNA"]]
DefaultAssay(lc5b)="RNA"










lc5b[["percent.mt"]] <- PercentageFeatureSet(lc5b, pattern = "^MT-")
VlnPlot(lc5b,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
multiplets <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
curva <- data.frame(multiplets,ncells)
ggplot(curva, aes(multiplets, ncells)) +
  geom_point() +
  geom_smooth(method = "lm")

fit <- lm(multiplets ~ ncells, curva)

model <- function(x){
  0.0007589*x + 0.0527214
}


perc <- model(ncol(lc5b))
q = (100 - perc)/100
feature_limit <- quantile(lc5b$nFeature_RNA, q)

lc5b <- subset(lc5b, subset = nFeature_RNA > 200 & nFeature_RNA < feature_limit & percent.mt < 10)

lc5b$patient="LC5"


lc5b$diagnosis="LongCovid"

lc5b$region="basal"
saveRDS(lc5b, file="LC5basal.rds")

lcvelo=merge(lc1a, y=c(lc1b, lc2a, lc2b, lc3a, lc3b, lc4, lc5a, lc5b))

saveRDS(lcvelo, file="lcvelo.rds")

lcvelocount=GetAssayData(lcvelo,slot="counts" )


lcvelo=CreateSeuratObject(lcvelocount)
lcvelo=UpdateGenesSeurat(lcvelo)

lcvelo=GetAssayData(lcvelo, slot="counts")
lcvelo=CreateSeuratObject(lcvelo, project = "lcvelo")

DotPlot(cedarsep, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(lcvelo, file="lcvelo.rds")


####### For data sets which are included GSE146981 were QC'ed the same as above, 
####so we will not show the QC step for these data sets, but just update the HGNC symbol.

load(("/Volumes/singlecell1/ipfall/ipfall2/ipf1883fep2.rds"))
meta=ipf1883fep@meta.data
count=GetAssayData(ipf1883fep, slot="counts")

a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ipf1883fep", colnames(count), sep="_")
a=colnames(count)


ipf1883fep=CreateSeuratObject(counts=count)


ipf1883fep=UpdateGenesSeurat(ipf1883fep)

ipf1883fep=GetAssayData(ipf1883fep, slot="counts")


ipf1883fep=CreateSeuratObject(counts=ipf1883fep, project ="ipf1883fep")



DotPlot(ipf1883fep, features=c("RACK1", "SELENOP", "SEPP1"))



saveRDS(ipf1883fep, file="ipf1883fep.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/ipf1880fep2.rds"))
VlnPlot(object=ipf1880fep,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

ipf1880fepmeta=ipf1880fep@meta.data
count=GetAssayData(ipf1880fep, slot="counts")

a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ipf1880fep", colnames(count), sep="_")
a=colnames(count)
ipf1880fep=CreateSeuratObject(counts=count)


ipf1880fep=UpdateGenesSeurat(ipf1880fep)

ipf1880fep=GetAssayData(ipf1880fep, slot="counts")
ipf1880fep=CreateSeuratObject(ipf1880fep, project = "ipf1880fep")

ipf1880fep@meta.data=ipf1880fepmeta

DotPlot(ipf1880fep, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf1880fep, file="ipf1880fep.rds")

load(("/Volumes/singlecell1/ipfall/ipfall2/ipf01bep18_2.rds"))
meta=ipf01bep18@meta.data
ipf01bep18count=GetAssayData(ipf01bep18, slot="counts")
VlnPlot(object=ipf01bep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf01bep18count)
colnames(ipf01bep18count)=str_replace(colnames(ipf01bep18count),  "-1", "")
a=colnames(ipf01bep18count)
colnames(ipf01bep18count)=paste("ipf01bep18", colnames(ipf01bep18count), sep="_")
a=colnames(ipf01bep18count)





ipf01bep18=CreateSeuratObject(counts=ipf01bep18count)
ipf01bep18@meta.data=meta

ipf01bep18=UpdateGenesSeurat(ipf01bep18)


ipf01bep18=GetAssayData(ipf01bep18, slot="counts")
ipf01bep18=CreateSeuratObject(ipf01bep18, project = "ipf01bep18")
DotPlot(ipf01bep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf01bep18, file="ipf01bep18.rds")




load(("/Volumes/singlecell1/ipfall/ipfall2/ipf02fep18_2.rds"))
meta=ipf02fep18@meta.data
count=GetAssayData(ipf02fep18, slot="counts")
VlnPlot(object=ipf02fep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ipf02fep18", colnames(count), sep="_")
a=colnames(count)





ipf02fep18=CreateSeuratObject(counts=count)
ipf02fep18@meta.data=meta

ipf02fep18=UpdateGenesSeurat(ipf02fep18)

ipf02fep18=GetAssayData(ipf02fep18, slot="counts")


ipf02fep18=CreateSeuratObject(counts=ipf02fep18, project ="ipf02fep18")





DotPlot(ipf02fep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf02fep18, file="ipf02fep18.rds")





load(("/Volumes/singlecell1/ipfall/ipfall2/ipf07fep17_2.rds"))
meta=ipf07fep17@meta.data
count=GetAssayData(ipf07fep17, slot="counts")
VlnPlot(object=ipf07fep17,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ipf07fep17", colnames(count), sep="_")
a=colnames(count)





ipf07fep17=CreateSeuratObject(counts=count)
ipf07fep17@meta.data=meta

ipf07fep17=UpdateGenesSeurat(ipf07fep17)

ipf07fep17=GetAssayData(ipf07fep17, slot="counts")


ipf07fep17=CreateSeuratObject(counts=ipf07fep17, project ="ipf07fep17")

DotPlot(ipf07fep17, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf07fep17, file="ipf07fep17.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/ca02nep_2.rds"))
meta=ca02nep@meta.data
count=GetAssayData(ca02nep, slot="counts")
VlnPlot(object=ca02nep,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ca02", colnames(count), sep="_")
a=colnames(count)





ca02nep=CreateSeuratObject(counts=count)
ca02nep@meta.data=meta

ca02nep=UpdateGenesSeurat(ca02nep)

ca02nep=GetAssayData(ca02nep, slot="counts")
ca02nep=CreateSeuratObject(ca02nep, project = "ca02nep")


DotPlot(ca02nep, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ca02nep, file="ca02nep.rds")


load(("/Volumes/singlecell1/ipfall/ipfall2/ca07nep_2.rds"))
meta=ca07nep@meta.data
count=GetAssayData(ca07nep, slot="counts")
VlnPlot(object=ca07nep,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ca07", colnames(count), sep="_")
a=colnames(count)





ca07nep=CreateSeuratObject(counts=count)
ca07nep@meta.data=meta

ca07nep=UpdateGenesSeurat(ca07nep)

ca07nep=GetAssayData(ca07nep, slot="counts")
ca07nep=CreateSeuratObject(ca07nep, project = "ca07nep")

DotPlot(ca07nep, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ca07nep, file="ca07nep.rds")
oad(("/Volumes/singlecell1/ipfall/ipfall2/ca12nep_2.rds"))
meta=ca12nep@meta.data
count=GetAssayData(ca12nep, slot="counts")
VlnPlot(object=ca12nep,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ca12nep", colnames(count), sep="_")
a=colnames(count)





ca12nep=CreateSeuratObject(counts=count)
ca12nep@meta.data=meta

ca12nep=UpdateGenesSeurat(ca12nep)

ca12nep=GetAssayData(ca12nep, slot="counts")
ca12nep=CreateSeuratObject(ca12nep, project = "ca12nep")




DotPlot(ca12nep, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ca12nep, file="ca12nep.rds")


load(("/Volumes/singlecell1/ipfall/ipfall2/ca16nep_2.rds"))
meta=ca16nep@meta.data
count=GetAssayData(ca16nep, slot="counts")
VlnPlot(object=ca16nep,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(count)
colnames(count)=str_replace(colnames(count),  "-1", "")
a=colnames(count)
colnames(count)=paste("ca16nep", colnames(count), sep="_")
a=colnames(count)





ca16nep=CreateSeuratObject(counts=count)
ca16nep@meta.data=meta

ca16nep=UpdateGenesSeurat(ca16nep)

ca16nep=GetAssayData(ca16nep, slot="counts")
ca16nep=CreateSeuratObject(ca16nep, project = "ca16nep")

DotPlot(ca16nep, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ca16nep, file="ca16nep.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/cc04dep17_2.rds"))
meta=cc04dep17@meta.data
cc04dep17count=GetAssayData(cc04dep17, slot="counts")
VlnPlot(object=cc04dep17,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(cc04dep17count)
colnames(cc04dep17count)=str_replace(colnames(cc04dep17count),  "-1", "")
a=colnames(cc04dep17count)
colnames(cc04dep17count)=paste("cc04dep17", colnames(cc04dep17count), sep="_")
a=colnames(cc04dep17count)





cc04dep17=CreateSeuratObject(counts=cc04dep17count)
cc04dep17@meta.data=meta

cc04dep17=UpdateGenesSeurat(cc04dep17)

cc04dep17=GetAssayData(cc04dep17, slot="counts")
cc04dep17=CreateSeuratObject(cc04dep17, project = "cc04dep17")

DotPlot(cc04dep17, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(cc04dep17, file="cc04dep17.rds")



oad(("/Volumes/singlecell1/ipfall/ipfall2/dd09dep18_2.rds"))
meta=dd09dep18@meta.data
dd09dep18count=GetAssayData(dd09dep18, slot="counts")
VlnPlot(object=dd09dep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(dd09dep18count)
colnames(dd09dep18count)=str_replace(colnames(dd09dep18count),  "-1", "")
a=colnames(dd09dep18count)
colnames(dd09dep18count)=paste("dd09dep18", colnames(dd09dep18count), sep="_")
a=colnames(dd09dep18count)





dd09dep18=CreateSeuratObject(counts=dd09dep18count)
dd09dep18@meta.data=meta

dd09dep18=UpdateGenesSeurat(dd09dep18)

dd09dep18=GetAssayData(dd09dep18, slot="counts")
dd09dep18=CreateSeuratObject(dd09dep18, project = "dd09dep18")

DotPlot(dd09dep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(dd09dep18, file="dd09dep18.rds")






load(("/Volumes/singlecell1/ipfall/ipfall2/dd10dep18_2.rds"))
meta=dd10dep18@meta.data
dd10dep18count=GetAssayData(dd10dep18, slot="counts")
VlnPlot(object=dd10dep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(dd10dep18count)
colnames(dd10dep18count)=str_replace(colnames(dd10dep18count),  "-1", "")
a=colnames(dd10dep18count)
colnames(dd10dep18count)=paste("dd10dep18", colnames(dd10dep18count), sep="_")
a=colnames(dd10dep18count)





dd10dep18=CreateSeuratObject(counts=dd10dep18count)
dd10dep18@meta.data=meta

dd10dep18=UpdateGenesSeurat(dd10dep18)


dd10dep18=GetAssayData(dd10dep18, slot="counts")
dd10dep18=CreateSeuratObject(dd10dep18, project = "dd10dep18")


DotPlot(dd10dep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(dd10dep18, file="dd10dep18.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/dd39dep17_2.rds"))
meta=dd39dep17@meta.data
dd39dep17count=GetAssayData(dd39dep17, slot="counts")
VlnPlot(object=dd39dep17,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(dd39dep17count)
colnames(dd39dep17count)=str_replace(colnames(dd39dep17count),  "-1", "")
a=colnames(dd39dep17count)
colnames(dd39dep17count)=paste("dd39dep17", colnames(dd39dep17count), sep="_")
a=colnames(dd39dep17count)





dd39dep17=CreateSeuratObject(counts=dd39dep17count)
dd39dep17@meta.data=meta

dd39dep17=UpdateGenesSeurat(dd39dep17)

dd39dep17=GetAssayData(dd39dep17, slot="counts")
dd39dep17=CreateSeuratObject(dd39dep17, project = "dd39dep17")

DotPlot(dd39dep17, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(dd39dep17, file="dd39dep17.rds")




load(("/Volumes/singlecell1/ipfall/ipfall2/dd49dep17_2.rds"))
meta=dd49dep17@meta.data
dd49dep17count=GetAssayData(dd49dep17, slot="counts")
VlnPlot(object=dd49dep17,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(dd49dep17count)
colnames(dd49dep17count)=str_replace(colnames(dd49dep17count),  "-1", "")
a=colnames(dd49dep17count)
colnames(dd49dep17count)=paste("dd49dep17", colnames(dd49dep17count), sep="_")
a=colnames(dd49dep17count)





dd49dep17=CreateSeuratObject(counts=dd49dep17count)
dd49dep17@meta.data=meta

dd49dep17=UpdateGenesSeurat(dd49dep17)
dd49dep17=GetAssayData(dd49dep17, slot="counts")
dd49dep17=CreateSeuratObject(dd49dep17, project = "dd49dep17")

DotPlot(dd49dep17, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(dd49dep17, file="dd49dep17.rds")













load(("/Volumes/singlecell1/ipfall/ipfall2/ipf12fep18_2.rds"))

ipf12fep18count=GetAssayData(ipf12fep18, slot="counts")
VlnPlot(object=ipf12fep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf12fep18count)
colnames(ipf12fep18count)=str_replace(colnames(ipf12fep18count),  "-1", "")
a=colnames(ipf12fep18count)
colnames(ipf12fep18count)=paste("ipf12fep18", colnames(ipf12fep18count), sep="_")
a=colnames(ipf12fep18count)





ipf12fep18=CreateSeuratObject(counts=ipf12fep18count)
DotPlot(ipf12fep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf12fep18, file="ipf12fep18.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/ipf07fep18_2.rds"))

ipf07fep18count=GetAssayData(ipf07fep18, slot="counts")
VlnPlot(object=ipf07fep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf07fep18count)
colnames(ipf07fep18count)=str_replace(colnames(ipf07fep18count),  "-1", "")
a=colnames(ipf07fep18count)
colnames(ipf07fep18count)=paste("ipf07fep18", colnames(ipf07fep18count), sep="_")
a=colnames(ipf07fep18count)





ipf07fep18=CreateSeuratObject(counts=ipf07fep18count)
DotPlot(ipf07fep18, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf07fep18, file="ipf07fep18.rds")




load(("/Volumes/singlecell1/ipfall/ipfall2/ipf05_19_2.rds"))

ipf05_19count=GetAssayData(ipf05_19, slot="counts")
VlnPlot(object=ipf05_19,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf05_19count)
colnames(ipf05_19count)=str_replace(colnames(ipf05_19count),  "-1_1", "")
a=colnames(ipf05_19count)
colnames(ipf05_19count)=paste("ipf05_19", colnames(ipf05_19count), sep="_")
a=colnames(ipf05_19count)





ipf05_19=CreateSeuratObject(counts=ipf05_19count)
DotPlot(ipf05_19, features=c("RACK1", "SELENOP", "SEPP1"))

saveRDS(ipf05_19, file="ipf05_19.rds")


load(("/Volumes/singlecell1/ipfall/ipfall2/ipf04fep18_2.rds"))

ipf04fep18count=GetAssayData(ipf04fep18, slot="counts")
VlnPlot(object=ipf04fep18,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf04fep18count)
colnames(ipf04fep18count)=str_replace(colnames(ipf04fep18count),  "-1", "")
a=colnames(ipf04fep18count)
colnames(ipf04fep18count)=paste("ipf04fep18", colnames(ipf04fep18count), sep="_")
a=colnames(ipf04fep18)


ipf04fep18=GetAssayData(ipf04fep18, slot="counts")
ipf04fep18=CreateSeuratObject(counts=ipf04fep18)


ipf04fep18=UpdateGenesSeurat(ipf04fep18)
ipf04fep18=GetAssayData(ipf04fep18, slot="counts")
ipf04fep18=CreateSeuratObject(ipf04fep18, project = "ipf04fep18")

ipf04fep18=CreateSeuratObject(counts=ipf04fep18)

DotPlot(ipf04fep18, features=c("RACK1", "SELENOP", "H3F3A", "H3-3A"))


saveRDS(ipf04fep18, file="ipf04fep18.rds")



load(("/Volumes/singlecell1/ipfall/ipfall2/ipf04_19_2.rds"))

DefaultAssay(ipf04_19)="RNA"
ipf04_19=GetAssayData(ipf04_19, slot="counts")
VlnPlot(object=ipf04_19,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf04_19)
colnames(ipf04_19)=str_replace(colnames(ipf04_19),  "-1", "")
a=colnames(ipf04_19)
colnames(ipf04_19)=paste("ipf04_19", colnames(ipf04_19), sep="_")
a=colnames(ipf04_19)





ipf04_19=CreateSeuratObject(counts=ipf04_19)



##QC check
VlnPlot(ipf04_19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


ipf04_19$patient="IPF04_19"

ipf04_19$diagnosis="IPF"

ipf04_19$hist="cystic"

ipf04_19$region="fibrotic"

save(ipf04_19, file="/Volumes/singlecell1/ipfall/ipfall2/ipf04_19_2.rds")




DotPlot(ipf04_19, features=c("RACK1", "SELENOP", "SEPP1"))


saveRDS(ipf04_19, file="ipf04_19.rds")







load(("/Volumes/singlecell1/ipfall/ipfall2/ipf03fep19_2.rds"))

ipf03fep19count=GetAssayData(ipf03fep19, slot="counts")
VlnPlot(object=ipf03fep19,c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)



a=colnames(ipf03fep19count)
colnames(ipf03fep19count)=str_replace(colnames(ipf03fep19count),  "-1", "")
a=colnames(ipf03fep19count)
colnames(ipf03fep19count)=paste("ipf03fep19", colnames(ipf03fep19count), sep="_")
a=colnames(ipf03fep19count)





ipf03fep19=CreateSeuratObject(counts=ipf03fep19count)


saveRDS(ipf03fep19, file="ipf03fep19.rds")




DotPlot(ipf03fep19, features=c("RACK1", "SELENOP", "SEPP1"))






DotPlot(cedarsep, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))


ipf05_19 =GetAssayData(ipf05_19 , slot="counts")

a=colnames(ipf05_19)

ipf05_19 =CreateSeuratObject(ipf05_19)
ipf05_19 =UpdateGenesSeurat(ipf05_19 )

ipf05_19 =GetAssayData(ipf05_19 , slot="counts")
ipf05_19 =CreateSeuratObject(ipf05_19 , project = "ipf05_19 ")

DotPlot(ipf05_19 , features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(ipf05_19 , file="ipf05_19 .rds")

ipf12fep18=GetAssayData(ipf12fep18, slot="counts")

a=colnames(ipf12fep18)

ipf12fep18=CreateSeuratObject(ipf12fep18)
ipf12fep18=UpdateGenesSeurat(ipf12fep18)

ipf12fep18=GetAssayData(ipf12fep18, slot="counts")
ipf12fep18=CreateSeuratObject(ipf12fep18, project = "ipf12fep18")

DotPlot(ipf12fep18, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(ipf12fep18, file="ipf12fep18.rds")

epmeta=ep@meta.data



ipf03fep19=GetAssayData(ipf03fep19, slot="counts")

a=colnames(ipf03fep19)

ipf03fep19=CreateSeuratObject(ipf03fep19)
ipf03fep19=UpdateGenesSeurat(ipf03fep19)

ipf03fep19=GetAssayData(ipf03fep19, slot="counts")
ipf03fep19=CreateSeuratObject(ipf03fep19, project = "ipf03fep19")

DotPlot(ipf03fep19, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(ipf03fep19, file="ipf03fep19.rds")



ipf07fep18=GetAssayData(ipf07fep18, slot="counts")

a=colnames(ipf07fep18)

ipf07fep18=CreateSeuratObject(ipf07fep18)
ipf07fep18=UpdateGenesSeurat(ipf07fep18)

ipf07fep18=GetAssayData(ipf07fep18, slot="counts")
ipf07fep18=CreateSeuratObject(ipf07fep18, project = "ipf07fep18")

DotPlot(ipf07fep18, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(ipf07fep18, file="ipf07fep18.rds")




ipf04_19=GetAssayData(ipf04_19, slot="counts")

a=colnames(ipf04_19)

ipf04_19=CreateSeuratObject(ipf04_19)
ipf04_19=UpdateGenesSeurat(ipf04_19)

ipf04_19=GetAssayData(ipf04_19, slot="counts")
ipf04_19=CreateSeuratObject(ipf04_19, project = "ipf04_19")

DotPlot(ipf04_19, features=c("RACK1", "SELENOP", "SEPP1", "H3F3A", "H3-3A"))

saveRDS(ipf04_19, file="ipf04_19.rds")




counts <- readMM("matrix.mtx")

# Read in `genes.tsv`
genes <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cells <- read_tsv("GSE135893_barcodes.tsv.gz", col_names = FALSE)
cell_ids=cells$X1
counts <- as(counts, "dgCMatrix")
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

meta=read.csv("metadata.csv", header=T, row.names = 1)

vand=CreateSeuratObject(counts = counts)
vand=UpdateGenesSeurat(vand)

vand@meta.data=meta
DotPlot(vand, features="RACK1", group.by="orig.ident")
DotPlot(vand, features="GNB2L1")
saveRDS(vand, file="vand.rds")
Idents(vand)="Diagnosis"
levels(vand)
vand2=subset(vand, idents=c("IPF", "Control" ))
all=merge(x=ipf1880fep, y=c(ipf1883fep, ipf01bep18, ipf02fep18, ipf07fep17, ipf03fep19, ipf04, ipf04fep18,ipf05,ipf07fep17 , ca02nep, ca07nep, ca12nep, ca16nep,
                               dd09dep18, cc04dep17, dd10dep18, dd39dep17, dd49dep17, ipf07fep18, ipf12fep18, lcvelo, vand2))

#### To remove low quality cells, we will only keep gene number greater than 200, percent.mt smaller than 10.
all=subset(all, subset = nFeature_RNA > 200 &  percent.mt < 10)

saveRDS(all, file="all.rds")













































