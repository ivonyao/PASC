library(Seurat)
library(ggplot2)
ep <- readRDS("~/ep.rds")
atm=read.csv("ATM signaling.csv", header=F)
atmgene=list(atm$V1)
ep=AddModuleScore(ep, features=atmgene, nbin=24, ctrl=100, name="ATM Signaling")

ep@meta.data[["ATM Signaling"]]=ep@meta.data[["ATM.Signaling1"]]
FeaturePlot(ep, features="ATM Signaling", split.by="diagnosis", min.cutoff = "q10", order=T)

eif2=read.csv("EIF2 list.csv", header=F)
eif2gene=list(eif2$V1)
ep=AddModuleScore(ep, features=eif2gene, nbin=24, ctrl=100, name="EIF2 Signaling")

ep@meta.data[["EIF2 Signaling"]]=ep@meta.data[["EIF2.Signaling1"]]


p53=read.csv("p53 signaling.csv", header=F)
p53gene=list(p53$V1)
ep=AddModuleScore(ep, features=p53gene, nbin=24, ctrl=100, name="p53 Signaling")

ep@meta.data[["p53 Signaling"]]=ep@meta.data[["p53.Signaling1"]]
FeaturePlot(ep, features="p53 Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")

Sirtuin=read.csv("Sirtuin.csv", header=F)
Sirtuingene=list(Sirtuin$V1)
ep=AddModuleScore(ep, features=Sirtuingene, nbin=24, ctrl=100, name="Sirtuin Signaling")

ep@meta.data[["Sirtuin Signaling"]]=ep@meta.data[["Sirtuin.Signaling1"]]
FeaturePlot(ep, features="Sirtuin Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")

saveRDS(ep, file="ep_modulescore.rds")


mTOR=read.csv("mTOR.csv", header=F)
mTORgene=list(mTOR$V1)
ep=AddModuleScore(ep, features=mTORgene, nbin=24, ctrl=100, name="mTOR Signaling")

ep@meta.data[["mTOR Signaling"]]=ep@meta.data[["mTOR.Signaling1"]]
FeaturePlot(ep, features="mTOR Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")
granzyma=read.csv("Granzyme A.csv", header=F)
granzymagene=list(granzyma$V1)
ep=AddModuleScore(ep, features=granzymagene, nbin=24, ctrl=100, name="Granzyme A Signaling")

ep@meta.data[["Granzyme A Signaling"]]=ep@meta.data[["Granzyme.A.Signaling1"]]
FeaturePlot(ep, features="Granzyme A Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")




ILK=read.csv("ILK.csv", header=F)
ILKgene=list(ILK$V1)
ep=AddModuleScore(ep, features=ILKgene, nbin=24, ctrl=100, name="ILK Signaling")

ep@meta.data[["ILK Signaling"]]=ep@meta.data[["ILK.Signaling1"]]
FeaturePlot(ep, features="ILK Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")

UPR=read.csv("UPR.csv", header=F)
UPRgene=list(UPR$V1)
ep=AddModuleScore(ep, features=UPRgene, nbin=24, ctrl=100, name="Unfolded protein response")

ep@meta.data[["Unfolded protein response"]]=ep@meta.data[["Unfolded.protein.response1"]]
FeaturePlot(ep, features="Unfolded protein response", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")


pru=read.csv("protein ubiquitination.csv", header=F)
prugene=list(pru$V1)
ep=AddModuleScore(ep, features=prugene, nbin=24, ctrl=100, name="Protein Ubiquitination Pathway")

ep@meta.data[["Protein Ubiquitination Pathway"]]=ep@meta.data[["Protein.Ubiquitination.Pathway1"]]
FeaturePlot(ep, features="Protein Ubiquitination Pathway", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")
Tgfbsig=read.csv("Tgfbsig.csv", header=F)
Tgfbsiggene=list(Tgfbsig$V1)
ep=AddModuleScore(ep, features=Tgfbsiggene, nbin=24, ctrl=100, name="TGF-β Signaling")

ep@meta.data[["TGF-β Signaling"]]=ep@meta.data[["TGF.β.Signaling1"]]



pf=read.csv("pulmonary fiborsis.csv", header=F)
pfgene=list(pf$V1)
ep=AddModuleScore(ep, features=pfgene, nbin=24, ctrl=100, name="Pulmonary Fibrosis Idiopathic Signaling Pathway")

ep@meta.data[["Pulmonary Fibrosis Idiopathic Signaling Pathway"]]=ep@meta.data[["Pulmonary.Fibrosis.Idiopathic.Signaling.Pathway1"]]
FeaturePlot(ep, features="Pulmonary Fibrosis Idiopathic Signaling Pathway", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")













FeaturePlot(ep, features="TGF-β Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")








fe=read.csv("Ferroptosis Signaling Pathway.csv", header=F)
fegene=list(fe$V1)
ep=AddModuleScore(ep, features=fegene, nbin=24, ctrl=100, name="Ferroptosis Signaling Pathway")

ep@meta.data[["Ferroptosis Signaling Pathway"]]=ep@meta.data[["Ferroptosis.Signaling.Pathway1"]]
FeaturePlot(ep, features="Ferroptosis Signaling Pathway", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



hif1a=read.csv("HIF1a signaling.csv", header=F)
hif1agene=list(hif1a$V1)
ep=AddModuleScore(ep, features=hif1agene, nbin=24, ctrl=100, name="HIF1α Signaling")

ep@meta.data[["HIF1α Signaling"]]=ep@meta.data[["HIF1α.Signaling1"]]

FeaturePlot(ep, features="HIF1α Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



telomere=read.csv("telomere.csv", header=F)
telomeregene=list(telomere$V1)
ep=AddModuleScore(ep, features=telomeregene, nbin=24, ctrl=100, name="Telomerase Signaling")

ep@meta.data[["Telomerase Signaling"]]=ep@meta.data[["Telomerase.Signaling1"]]
FeaturePlot(ep, features="Telomerase Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



Cholesterol=read.csv("cholestorl super pathway.csv", header=F)
Cholesterolgene=list(Cholesterol$V1)
ep=AddModuleScore(ep, features=Cholesterolgene, nbin=24, ctrl=100, name="Superpathway.of.Cholesterol.Biosynthesis1")

ep@meta.data[["Superpathway of Cholesterol Biosynthesis"]]=ep@meta.data[["Superpathway.of.Cholesterol.Biosynthesis1"]]
FeaturePlot(ep, features="Superpathway of Cholesterol Biosynthesis", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")




mt=read.csv("Mitochondria_dysfuntion.csv", header=F)
mtgene=list(mt$V1)
ep=AddModuleScore(ep, features=mtgene, nbin=24, ctrl=100, name="Mitochondrial Dysfunction")

ep@meta.data[["Mitochondrial Dysfunction"]]=ep@meta.data[["Mitochondrial.Dysfunction1"]]
FeaturePlot(ep, features="Mitochondrial Dysfunction", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



se=read.csv("senescence.csv", header=F)
segene=list(se$V1)
ep=AddModuleScore(ep, features=segene, nbin=24, ctrl=100, name="Senescence Pathway")

ep@meta.data[["Senescence Pathway"]]=ep@meta.data[["Senescence.Pathway1"]]
FeaturePlot(ep, features="Senescence Pathway", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



levels(ep)


DimPlot(ep)





HER=read.csv("Her-2.csv", header=F)
HERgene=list(HER$V1)
ep=AddModuleScore(ep, features=HERgene, nbin=24, ctrl=100, name="HER-2 Signaling in Breast Cancer")

ep@meta.data[["HER-2 Signaling in Breast Cancer"]]=ep@meta.data[["HER.2.Signaling.in.Breast.Cancer1"]]
FeaturePlot(ep, features="HER-2 Signaling in Breast Cancer", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")



nrf2=read.csv("NRF2-mediated Oxidative Stress Response.csv", header=F)
nrf2gene=list(nrf2$V1)
ep=AddModuleScore(ep, features=nrf2gene, nbin=24, ctrl=100, name="NRF2-mediated Oxidative Stress Response")

ep@meta.data[["NRF2-mediated Oxidative Stress Response"]]=ep@meta.data[["NRF2.mediated.Oxidative.Stress.Response1"]]
FeaturePlot(ep, features="NRF2-mediated Oxidative Stress Response", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")


rho=read.csv("RHOA.csv", header=F)
rhogene=list(rho$V1)
ep=AddModuleScore(ep, features=rhogene, nbin=24, ctrl=100, name="RHOA Signaling")

ep@meta.data[["RHOA Signaling"]]=ep@meta.data[["RHOA.Signaling1"]]
FeaturePlot(ep, features="RHOA Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")


hippo=read.csv("Hippo sig.csv", header=F)
hippogene=list(hippo$V1)
ep=AddModuleScore(ep, features=hippogene, nbin=24, ctrl=100, name="HIPPO Signaling")

ep@meta.data[["HIPPO Signaling"]]=ep@meta.data[["HIPPO.Signaling1"]]
FeaturePlot(ep, features="HIPPO Signaling", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")


saveRDS(ep, file="ep_modulescore.rds")
levels(ep)

cs=read.csv("csgene_human2.csv", header=F)
csgene=list(cs$V1)
ep=AddModuleScore(ep, features=csgene, nbin=24, ctrl=100, name="Cellular Senescence")

ep@meta.data[["Cellular Senescence"]]=ep@meta.data[["Cellular.Senescence1"]]
FeaturePlot(ep, features="Cellular Senescence", split.by="diagnosis", min.cutoff = "q10", order=T, max.cutoff = "q95")

ep=RenameIdents(ep, "KRT17+ KRT15+ Basal" ="Basal", "KRT17+ KRT15- Basal"="Basal" )

levels(ep)=c("AT1" ,     "AT2" , "Basal" ,    "Ciliated" ,  "Club",               "Cycling Epithelial", "Goblet" ,            "SCGB3A2+",           "Serous"  ,           "Transitional"   )

ep$celltype5=Idents(ep)

ep$cell_dia5=paste(ep$celltype5, ep$diagnosis, sep="-")


Idents(ep)=ep$cell_dia5
ave=AverageExpression(ep, return.seurat =T)


meta=ep@meta.data

write.csv(meta, "meta_score.csv")
matrix=FetchData(ep, vars=c("ATM Signaling",	"EIF2 Signaling",	"p53 Signaling",	"Sirtuin Signaling",	"mTOR Signaling",	"Granzyme A Signaling",	"ILK Signaling",	"Unfolded protein response",	"Protein Ubiquitination Pathway",	"TGF-β Signaling", 	"Pulmonary Fibrosis Idiopathic Signaling Pathway", 	"Ferroptosis Signaling Pathway", 	"HIF1α Signaling", 	"Telomerase Signaling",		"Superpathway of Cholesterol Biosynthesis",	"Mitochondrial Dysfunction",	"Cellular Senescence", 	"RHOA Signaling", 	"NRF2-mediated Oxidative Stress Response", 	"HER-2 Signaling in Breast Cancer",		"HIPPO Signaling"))


mat=t(matrix)

mats=CreateSeuratObject(mat)

mats$celldia=ep$celltype_diag

mats$celldia5=ep$cell_dia5
levels(mats)
Idents(mats)=mats$celldia
avemat=AverageExpression(mats, return.seurat = T)

avematrix=GetAssayData(avemat, slot="data")

am=t(avematrix)

am1=am[order(row.names(am)), order(colnames(am))]
am2=t(am1)  


Idents(mats)=mats$celldia5
avemat2=AverageExpression(mats, return.seurat = T)

amtrix=GetAssayData(avemat2, slot="data")

amtrix2=amtrix[ order(row.names(amtrix)), order(colnames(amtrix))]



library(pheatmap)
  at1=am2[ , 1:3]
  pheatmap(at1, scale="row", cluster_cols=F, cluster_rows = F)
  pheatmap(at1, scale="row", cluster_cols=F)
  
  
  pheatmap(at2, scale="row", cluster_cols=F , cluster_rows = F)
  
  
  pheatmap(at2, scale="row", cluster_cols=F )
  


am3=am2






vec=c("ATM Signaling", "Unfolded protein response", "HIPPO Signaling", 	"Telomerase Signaling",  "Pulmonary Fibrosis Idiopathic Signaling Pathway", "Protein Ubiquitination Pathway", 	"TGF-β Signaling", "p53 Signaling", "Cellular Senescence", "HIF1α Signaling",  "RHOA Signaling","Ferroptosis Signaling Pathway",  "ILK Signaling","EIF2 Signaling",		"Sirtuin Signaling",	"mTOR Signaling",	"Granzyme A Signaling",				 				"Superpathway of Cholesterol Biosynthesis",	"Mitochondrial Dysfunction",	 		"NRF2-mediated Oxidative Stress Response", 	"HER-2 Signaling in Breast Cancer")

am4=am3[order(match(rownames(am3), vec)), , drop = T]

at1=am4[ , 1:3]

at1 <- at1[, c(1, 3, 2)]



pheatmap(at1, scale="row", cluster_cols=F, cluster_rows = F)
pheatmap(at1, scale="row", cluster_cols=F)
pheatmap(at1, main = "AT1", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )

at2=am4[ , 4:6]
at2 <- at2[, c(1, 3, 2)]
pheatmap(at2, scale="row", cluster_cols=F , cluster_rows = F)
pheatmap(at2, main = "AT2", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )

pheatmap(at2, scale="row", cluster_cols=F )


cili=am4[ , 7:9]
cili <- cili[, c(1, 3, 2)]

pheatmap(cili, main = "Ciliated", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF" ))

club=am4[ , 10:12]

club <- club[, c(1, 3, 2)]
pheatmap(club, main = "Club", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )



gob=am4[ , 16:18]
gob <- gob[, c(1, 3, 2)]

pheatmap(gob, main = "Goblet", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )

k17=am4[ , 19:21]
k17 <- k17[, c(1, 3, 2)]
pheatmap(k17, main = "KRT17+ KRT15- Basal", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )


k15=am4[ , 22:24]
k15 <- k15[, c(1, 3, 2)]

pheatmap(k15, main = "KRT17+ KRT15+ Basal", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )



scgb=am4[ , 25:27]
scgb <- scgb[, c(1, 3, 2)]
pheatmap(scgb, main = "SCGB3A2+", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )


library(pheatmap)
tran=am4[ , 30:32]
tran <- tran[, c(1, 3, 2)]
pheatmap(tran, main = "Transitional", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )




am5=amtrix2[order(match(rownames(am3), vec)), , drop = T]


basal=am5[ , 7:9]
basal <- basal[, c(1, 3, 2)]
pheatmap(basal, main = "Basal", scale="row", cluster_cols=F , cluster_rows = F, angle_col=0,  labels_col=c("Control", "PASC", "IPF") )

    