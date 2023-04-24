library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(velocyto.R)
library(ggplot2)
library(igraph)
library(stringr)



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