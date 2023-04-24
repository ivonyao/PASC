

library(Seurat)

ep <- readRDS("~/ep.rds")
at2=subset(ep, idents="AT2")



Idents(at2)="diagnosis"


levels(at2)
at2_PASC_vs_control=FindMarkers(at2, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15,  test.use ="MAST")
saveRDS(at2_PASC_vs_control, file="at2_PASC_vs_control.rds")
write.csv(at2_PASC_vs_control, "at2_PASC_vs_control.csv")


at2_IPF_vs_control=FindMarkers(at2, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(at2_IPF_vs_control, file="at2_IPF_vs_control.rds")

write.csv(at2_IPF_vs_control, "at2_IPF_vs_control.csv")


at2_PASC_vs_IPF=FindMarkers(at2, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(at2_PASC_vs_IPF, file="at2_PASC_vs_IPF.rds")
write.csv(at2_PASC_vs_IPF, "at2_PASC_vs_IPF.csv")

saveRDS(at2, file="at2.rds")

at1=subset(ep, idents="AT1")





Idents(at1)="diagnosis"


levels(at1)
at1_PASC_vs_control=FindMarkers(at1, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(at1_PASC_vs_control, file="at1_PASC_vs_control.rds")
write.csv(at1_PASC_vs_control, "at1_PASC_vs_control.csv")


at1_IPF_vs_control=FindMarkers(at1, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(at1_IPF_vs_control, file="at1_IPF_vs_control.rds")

write.csv(at1_IPF_vs_control, "at1_IPF_vs_control.csv")


at1_PASC_vs_IPF=FindMarkers(at1, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(at1_PASC_vs_IPF, file="at1_PASC_vs_IPF.rds")
write.csv(at1_PASC_vs_IPF, "at1_PASC_vs_IPF.csv")


Ciliated=subset(ep, idents="Ciliated")

Idents(Ciliated)="diagnosis"


levels(Ciliated)
Ciliated_PASC_vs_control=FindMarkers(Ciliated, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Ciliated_PASC_vs_control, file="Ciliated_PASC_vs_control.rds")
write.csv(Ciliated_PASC_vs_control, "Ciliated_PASC_vs_control.csv")


Ciliated_IPF_vs_control=FindMarkers(Ciliated, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Ciliated_IPF_vs_control, file="Ciliated_IPF_vs_control.rds")

write.csv(Ciliated_IPF_vs_control, "Ciliated_IPF_vs_control.csv")


Ciliated_PASC_vs_IPF=FindMarkers(Ciliated, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Ciliated_PASC_vs_IPF, file="Ciliated_PASC_vs_IPF.rds")
write.csv(Ciliated_PASC_vs_IPF, "Ciliated_PASC_vs_IPF.csv")

saveRDS(Ciliated, file="Ciliated.rds")





basal=subset(ep, idents=c( "KRT17+ KRT15- Basal", "KRT17+ KRT15+ Basal"))




Idents(basal)="diagnosis"


levels(basal)
basal_PASC_vs_control=FindMarkers(basal, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(basal_PASC_vs_control, file="basal_PASC_vs_control.rds")
write.csv(basal_PASC_vs_control, "basal_PASC_vs_control.csv")


basal_IPF_vs_control=FindMarkers(basal, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(basal_IPF_vs_control, file="basal_IPF_vs_control.rds")

write.csv(basal_IPF_vs_control, "basal_IPF_vs_control.csv")


basal_PASC_vs_IPF=FindMarkers(basal, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(basal_PASC_vs_IPF, file="basal_PASC_vs_IPF.rds")
write.csv(basal_PASC_vs_IPF, "basal_PASC_vs_IPF.csv")


saveRDS(basal, file="basal.rds")



Goblet=subset(ep, idents=c( "Goblet"))


levels(ep)


Idents(Goblet)="diagnosis"


levels(Goblet)
Goblet_PASC_vs_control=FindMarkers(Goblet, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Goblet_PASC_vs_control, file="Goblet_PASC_vs_control.rds")
write.csv(Goblet_PASC_vs_control, "Goblet_PASC_vs_control.csv")


Goblet_IPF_vs_control=FindMarkers(Goblet, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Goblet_IPF_vs_control, file="Goblet_IPF_vs_control.rds")

write.csv(Goblet_IPF_vs_control, "Goblet_IPF_vs_control.csv")


Goblet_PASC_vs_IPF=FindMarkers(Goblet, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST")
saveRDS(Goblet_PASC_vs_IPF, file="Goblet_PASC_vs_IPF.rds")
write.csv(Goblet_PASC_vs_IPF, "Goblet_PASC_vs_IPF.csv")


saveRDS(Goblet, file="Goblet.rds")


SCGB3A2=subset(ep, idents=c( "SCGB3A2+"))
levels(ep)


Idents(SCGB3A2)="diagnosis"


levels(SCGB3A2)
SCGB3A2_PASC_vs_control=FindMarkers(SCGB3A2, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(SCGB3A2_PASC_vs_control, file="SCGB3A2_PASC_vs_control.rds")
write.csv(SCGB3A2_PASC_vs_control, "SCGB3A2_PASC_vs_control.csv")


SCGB3A2_IPF_vs_control=FindMarkers(SCGB3A2, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(SCGB3A2_IPF_vs_control, file="SCGB3A2_IPF_vs_control.rds")

write.csv(SCGB3A2_IPF_vs_control, "SCGB3A2_IPF_vs_control.csv")


SCGB3A2_PASC_vs_IPF=FindMarkers(SCGB3A2, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(SCGB3A2_PASC_vs_IPF, file="SCGB3A2_PASC_vs_IPF.rds")
write.csv(SCGB3A2_PASC_vs_IPF, "SCGB3A2_PASC_vs_IPF.csv")
saveRDS(SCGB3A2, file="SCGB3A2.rds")




Club=subset(ep, idents=c( "Club"))


levels(ep)


Idents(Club)="diagnosis"


levels(Club)
Club_PASC_vs_control=FindMarkers(Club, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Club_PASC_vs_control, file="Club_PASC_vs_control.rds")
write.csv(Club_PASC_vs_control, "Club_PASC_vs_control.csv")


Club_IPF_vs_control=FindMarkers(Club, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Club_IPF_vs_control, file="Club_IPF_vs_control.rds")

write.csv(Club_IPF_vs_control, "Club_IPF_vs_control.csv")


Club_PASC_vs_IPF=FindMarkers(Club, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Club_PASC_vs_IPF, file="Club_PASC_vs_IPF.rds")
write.csv(Club_PASC_vs_IPF, "Club_PASC_vs_IPF.csv")




Transitional=subset(ep, idents=c( "Transitional"))


levels(ep)


Idents(Transitional)="diagnosis"


levels(Transitional)
Transitional_PASC_vs_control=FindMarkers(Transitional, ident.1 = "PASC", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Transitional_PASC_vs_control, file="Transitional_PASC_vs_control.rds")
write.csv(Transitional_PASC_vs_control, "Transitional_PASC_vs_control.csv")


Transitional_IPF_vs_control=FindMarkers(Transitional, ident.1 = "IPF", ident.2="Control", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Transitional_IPF_vs_control, file="Transitional_IPF_vs_control.rds")

write.csv(Transitional_IPF_vs_control, "Transitional_IPF_vs_control.csv")


Transitional_PASC_vs_IPF=FindMarkers(Transitional, ident.1 = "PASC", ident.2="IPF", logfc.threshold = 0.15, test.use ="MAST", verbose=T)
saveRDS(Transitional_PASC_vs_IPF, file="Transitional_PASC_vs_IPF.rds")
write.csv(Transitional_PASC_vs_IPF, "Transitional_PASC_vs_IPF.csv")





























































































































































