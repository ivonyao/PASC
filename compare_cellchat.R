library(CellChat)
library(patchwork)
library(ggplot2)
library(viridis)
Control <- readRDS("/Volumes/T7/longcovid/Gene_symbol_update/refup_duprm/cellchat4/control/controlsub_cellchat0.23.rds")
IPF <- readRDS("/Volumes/T7/longcovid/Gene_symbol_update/refup_duprm/cellchat4/ipf/IPFsub_cellchat0.23.rds")
PASC <- readRDS("/Volumes/T7/longcovid/Gene_symbol_update/refup_duprm/cellchat4/pasc/PASCsub_cellchat0.23.rds")

Control <- updateCellChat(Control)
PASC <- updateCellChat(PASC)
IPF <- updateCellChat(IPF)
object.list <- list(Control = Control, PASC = PASC, IPF=IPF)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  555 signaling genes.
#>  7563 cells. 
#> CellChat analysis of single cell RNA-seq data!



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2








par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2




weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 4)


rankSimilarity(cellchat, type = "functional", comparison2 = c(1, 2, 3))

?rankSimilarity

gg1 <- rankNet(cellchat, mode = "comparison",  comparison = c(1, 2, 3), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison",  comparison = c(1, 2, 3), stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)

i = 1
pathway.union <- union(object.list[[i]]@netP[["pathways"]], object.list[[i+1]]@netP[["pathways"]])

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
ht3=netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6)

draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))




ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
ht3=netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6)

draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
ht3=netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 6, color.heatmap = "OrRd")

draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
levels(Control@idents)
netVisual_bubble(cellchat, sources.use = c(12, 16), targets.use = c(1:9),  comparison = c(1, 2), angle.x = 45)


a=netVisual_bubble(cellchat, sources.use = c(1, 2, 4, 5, 6, 7, 8, 9, 10), targets.use = c(3, 11),signaling = "TGFb",  comparison = c(1, 2, 3), angle.x = 45)
a=a[["data"]]
a$dataset=factor(a$dataset, levels=c("Control", "PASC", "IPF"))
ggplot(a, aes(x=dataset, y=group.names, size=pval   ))+geom_point(aes(colour = prob.original)) +
  scale_colour_gradient2(low = "blue",
                         mid = "green",
                         high = "red", na.value = "lightgrey")

a$prob2=a$prob*100
a=a[order(a$receptor), ]

c=a[17:70, ]

ggplot(c, aes(x=dataset, y=group.names, size=pval   ))+geom_point(aes(colour = prob2)) + scale_color_viridis(option="H")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

levels(cellchat@idents)

netVisual_bubble(cellchat, sources.use = c(11, 12, 13, 16, 18, 17, 19), targets.use = c(9),signaling = "TGFb",  comparison = c(1, 2, 3), angle.x = 45)



gg1 <- netVisual_bubble(cellchat, sources.use = 16, targets.use = c(1:9),  comparison = c(1, 2, 3), max.dataset = 1, title.name = "Increased signaling in Control", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use =16, targets.use = c(1:9),  comparison = c(1, 2, 3), max.dataset = 2, title.name = "Increased signaling in PASC", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg3<- netVisual_bubble(cellchat, sources.use =16, targets.use = c(1:9),  comparison = c(1, 2, 3), max.dataset = 3, title.name = "Increased signaling in IPF", angle.x = 45, remove.isolate = T)


gg1 + gg2+gg3

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PASC"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "PASC",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Control",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 16, targets.use = c(1:9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 16, targets.use = c(1:9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2




?netVisual_aggregate
pathways.show <- c("TGFb") 
pdf(file ="tgfbcho.pdf", width = 15, height =5)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show , signaling.name = paste(pathways.show, names(object.list)[i]), layout = "chord")
}
dev.off()

pdf(file ="tgfbLR.pdf", width = 15, height =5)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netAnalysis_signalingRole_network(object.list[[i]], signaling = pathways.show)
}
dev.off()




saveRDS(cellchat, file="comparesubset.rds")
saveRDS(object.list, file="object.listsubset.rds")
