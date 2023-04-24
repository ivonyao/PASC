install.packages('NMF')
devtools::install_github("jokergoo/circlize")
library(NMF)
devtools::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
devtools::install_github("sqjin/cellchat")
library(CellChat)
library(Seurat)

epim <- readRDS("/Volumes/T7/longcovid/Gene_symbol_update/refup_duprm/epim3.rds")

epim=subset(epim, idents="Serous", invert=T)
levels(epim)

Idents(epim)="diagnosis"
epim=RenameIdents(epim, "LongCovid"="PASC")
levels(epim)=c(   "Control", "PASC", "IPF")

epim$diagnosis=Idents(epim)




PASC=subset(epim, idents="PASC")



data=GetAssayData(PASC, assay = "RNA", slot = "data" )

Idents(PASC)="celltype5"

labels <- Idents(PASC)
table(PASC$diagnosis)

metadata=PASC@meta.data

meta <- data.frame(group = labels, row.names = names(labels))
unique(meta$group)

cellchat <- createCellChat(object = data, meta = meta, group.by = "group")

cellchat <- addMeta(cellchat, meta = metadata)
cellchat <- setIdent(cellchat, ident.use = "celltype5")

groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)





# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB



# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) 
options(future.globals.maxSize= 2000000000)

#> do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to PASC forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human)


?computeCommunProb



cellchat<- computeCommunProb(cellchat, type ="truncatedMean",  trim = 0.23)
#> triMean is used for caculating the average gene expression per cell group. 

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


df.net <- subsetCommunication(cellchat) 
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net2 <- subsetCommunication(cellchat, signaling = c( "TGFb")) 
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)





groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
levels(cellchat@idents)
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,10) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")







par(mfrow=c(1,1))

pdf(file ="tgfbcho1.pdf", width = 10, height =10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()



# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object



# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:22), remove.isolate = FALSE)
#> Comparing communications on a single object
#> 
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:23), remove.isolate = FALSE)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
saveRDS(cellchat, file="PASC_cellchat0.23.rds")





