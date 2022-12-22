setwd("/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/")
sham.data <- Read10X(data.dir = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/filtered_feature_bc_matrix_sham/")
sham <- CreateSeuratObject(counts = sham.data, project = "sham.meninges", min.cells = 3, min.features = 200)

TBI.data <- Read10X(data.dir = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/filtered_feature_bc_matrix_TBI/")
TBI <- CreateSeuratObject(counts = TBI.data, project = "TBI.meninges", min.cells = 3, min.features = 200)

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)

meninges.all<- merge(x = sham,
                     y = TBI,
                     add.cell.ids = c("Sham", "TBI"),
                     project = "meninges1wpi")
head(colnames(meninges.all))

meninges.all[["percent.mt"]] <- PercentageFeatureSet(meninges.all, pattern = "^mt-")
meninges.all[["percent.hemo"]] <- PercentageFeatureSet(meninges.all, pattern = "^Hbb-")

qc.featurerna.meninges.all <- VlnPlot(meninges.all, features = "nFeature_RNA") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")

qc.nCountrna.meninges.all <- VlnPlot(meninges.all, features = "nCount_RNA") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")

qc.percentMt.meninges.all <- VlnPlot(meninges.all, features = "percent.mt") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")

qc.percentHemo.meninges.all <- VlnPlot(meninges.all, features = "percent.hemo") + NoLegend() + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")

qc.count.featurerna.scatter.meninges.all <- FeatureScatter(meninges.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1)) + xlab("Total RNA reads") + 
  ylab("Number of expressed genes")

qc.count.percentmt.scatter.meninges.all <- FeatureScatter(meninges.all, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1)) + xlab("Total RNA reads") + 
  ylab("Percent Mitochondrial Genes")

qc.count.percenthemo.scatter.meninges.all <- FeatureScatter(meninges.all, feature1 = "nCount_RNA", feature2 = "percent.hemo")+
  theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1)) + xlab("Total RNA reads") + ylab("Percent Hemoglobin Beta") 


pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.QCgraphs.meninges.all.pdf")
(qc.featurerna.meninges.all | qc.nCountrna.meninges.all | qc.percentMt.meninges.all | qc.percentHemo.meninges.all) / (qc.count.featurerna.scatter.meninges.all | qc.count.percentmt.scatter.meninges.all | qc.count.percenthemo.scatter.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.featurena.meninges.all.pdf")
(qc.featurerna.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.nCountrna.meninges.all.pdf")
(qc.nCountrna.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.percentMt.meninges.all.pdf")
(qc.percentMt.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.percentHemo.meninges.all.pdf")
(qc.percentHemo.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.count.featurerna.scatter.meninges.all.pdf")
(qc.count.featurerna.scatter.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.count.percentmt.scatter.meninges.all.pdf")
(qc.count.percentmt.scatter.meninges.all)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.qc.count.percenthemo.scatter.meninges.all.pdf")
(qc.count.percenthemo.scatter.meninges.all)
dev.off()

meninges.all <- subset(meninges.all, subset = nFeature_RNA >= 150 & nFeature_RNA <= 5000 & 
                         percent.mt < 20 & percent.hemo < 5)

meninges.all <- NormalizeData(meninges.all, normalization.method = "LogNormalize", scale.factor = 10000)
meninges.all <- FindVariableFeatures(meninges.all, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.all), 10)
top10

plot1 <- VariableFeaturePlot(meninges.all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.variablefeatures.meninges.all.pdf")
(plot2)
dev.off()

all.genes.meninges.all <- rownames(meninges.all)
meninges.all <- ScaleData(meninges.all, features = all.genes.meninges.all, vars.to.regress = "percent.mt")

meninges.all <- RunPCA(meninges.all, features = VariableFeatures(object = meninges.all))
print(meninges.TBI[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.Vizdimreduction.pca.meninges.all.pdf")
VizDimLoadings(meninges.all, dims = 1:2, reduction = "pca")
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.dimplot.pca.meninges.all.pdf")
DimPlot(meninges.all, reduction = "pca")
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.Dimheatmap.10dim.meninges.all.pdf")
DimHeatmap(meninges.all, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.Dimheatmap.1dim.meninges.all.pdf")
DimHeatmap(meninges.all, dims = 1, cells = 500, balanced = TRUE)
dev.off()

meninges.all <- JackStraw(meninges.all, num.replicate = 100)
meninges.all <- ScoreJackStraw(meninges.all, dims = 1:20)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.jackstrawdiminsionality.all.pdf")
JackStrawPlot(meninges.all, dims = 1:20)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.01.2020.elbowplotdiminsionality.all.pdf")
ElbowPlot(meninges.all)
dev.off()

meninges.all <- FindNeighbors(meninges.all, dims = 1:15)
meninges.all <- FindClusters(meninges.all, resolution = 0.5)
head(Idents(meninges.all), 5)

meninges.all <- RunUMAP(meninges.all, dims = 1:15)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', cols = c("cornflowerblue", "darkseagreen2"), order = c("TBI", "Sham"))

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.UMAP.meninges.all.pdf", width = 10, height =8)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', cols = c("orchid1","darkturquoise"), order = c("Sham", "TBI"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.UMAP.meninges.all.split.pdf", width = 14, height =8)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = c("orchid1","darkturquoise"), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.UMAP.meninges.all.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.all, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.UMAP.meninges.all.cluster_colors_NOLABEL.pdf", width = 10, height =8)
DimPlot(meninges.all, reduction = "umap", label = FALSE, repel = FALSE)
dev.off()

View(meninges.all)
saveRDS(meninges.all, file = "meninges.all.122120.rds")

meninges.all.markers <- FindAllMarkers(meninges.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.gene <- (meninges.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))
View(top.gene)

saveRDS(meninges.all.markers, file = "meninges.all.markers.122120.rds")
write.table(meninges.all.markers, file = "meninges_all_markers_122120.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####Assign Cluster Names####
meninges.all.markers<- readRDS('meninges.all.markers.122120.rds')
meninges.all <- readRDS('meninges.all.122120.rds')

meninges.all[["OriginalClusters"]] <- Idents(object = meninges.all)
new.cluster.ids <- c("Fibroblasts","Activated Macrophages 1","Endothelial Cells 1","Activated Macrophages 2","CD3+ T Cells","B Cells 1","Dendritic Cells","Activated T Cells","B Cells 2","Endothelial Cells 2","Choroid Plexus","Immature/Differentiating B Cells","NK Cells","Pineal Gland Cells","Proliferating Cells","Pericytes","Neural Crest Cells","Clotting Related","Plasmacytoid Dendritic Cells","Macrophages 3","Neutrophils")
names(new.cluster.ids) <- levels(meninges.all)
meninges.all.names <- RenameIdents(meninges.all, new.cluster.ids)

my_levels <- c("Activated Macrophages 1","Activated Macrophages 2","Macrophages 3","B Cells 1","B Cells 2","Immature/Differentiating B Cells","CD3+ T Cells","Activated T Cells","NK Cells","Dendritic Cells","Plasmacytoid Dendritic Cells","Neutrophils","Proliferating Cells","Fibroblasts","Endothelial Cells 1","Endothelial Cells 2","Pericytes","Choroid Plexus","Pineal Gland Cells","Neural Crest Cells", "Clotting Related")
meninges.all.names@active.ident <- factor(x=meninges.all.names@active.ident, levels = my_levels)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.meninges.all.cluster_colors_withlabels.pdf", width = 14, height =8)
DimPlot(meninges.all.names, reduction = "umap")
dev.off()

saveRDS(meninges.all.names, file = "meninges.all.names.011220.rds")

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.11.2021.UMAP.meninges.all.split.withnames.pdf", width = 16, height =8)
DimPlot(meninges.all.names, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham", "TBI"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.meninges.all.names.WITHREPEL.pdf", width = 14, height =8)
DimPlot(meninges.all.names, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)
dev.off()

meninges.all.names[["OriginalClusters"]] <- Idents(object = meninges.all.names)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")
names(new.cluster.ids) <- levels(meninges.all.names)
meninges.all.numbers <- RenameIdents(meninges.all.names, new.cluster.ids)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.meninges.all.cluster_colors_withnumbers.pdf", width = 10, height =8)
DimPlot(meninges.all.numbers, reduction = "umap", label = TRUE)
dev.off()

####find the markers of each cluster####
meninges.all.markers.names <- FindAllMarkers(meninges.all.names, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(meninges.all.markers.names)
top.gene.names.all <- (meninges.all.markers.names %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))
View(top.gene.names.all)
write.table(top.gene.names.all, file = "top_gene_names_bycluster_all_080320.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.endothelial_cells.pdf", width = 10, height =8)
DimPlot(meninges.all.numbers, reduction = "UMAP", label = TRUE, idents = c("16","17"))
dev.off()

Cluster.fibroblasts.only <- subset(x = meninges.all.names, idents = c('Fibroblasts'))
pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.fibroblasts.pdf", width = 10, height =8)
DimPlot(Cluster.fibroblasts.only, reduction = "umap", label = TRUE)
dev.off()


meninges.endo <- FindNeighbors(Cluster.endothelial.only, dims = 1:15)
meninges.endo <- FindClusters(Cluster.endothelial.only, resolution = 0.5)
meninges.endo.markers <- FindAllMarkers(meninges.endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster.endothelial.only <- subset(x = meninges.all.names, idents = c('Endothelial Cells 1', 'Endothelial Cells 2'))
pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//12.21.2020.UMAP.endothelial.pdf", width = 10, height =8)
DimPlot(Cluster.endothelial.only, reduction = "tsne", label = TRUE)
dev.off()

main.feature.plots <- (FeaturePlot(meninges.all.names, 
                                reduction = "umap", 
                                pt.size = 3,
                                features = c("Nkr1")))
main.feature.plots

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.31.2021.Tacr1.feature.pdf", width = 8, height =8)
main.feature.plots
dev.off()

####make dot plot####

print(top.2.genes[,7])

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.11.2021.dotplot_topgenes_withlabels.all.pdf", width = 10, height = 6)
DotPlot(meninges.all.names,
        features = c("C1qb","Pf4","Fscn1","Ms4a1","Iglc2","Vpreb3","Cd3g","Il7r","Gzma","Cd209a","Cox6a2","S100a9","Mki67","Col1a1","Igfbp3","Vwf","Rgs5","Enpp2","Chgb","Mpz","Plac8"),
        assay = "RNA", 
        cols = c("grey80", "red"), col.min = 0) + FontSize(y.title = 12, x.text = 12, y.text =12) +RotatedAxis() + xlab("") +
  ylab("Cluster Number") + scale_color_gradientn(colors = rev(rainbow(64, start = 0, end =0.7)))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//08.03.2020.dotplot_topgenes_withlabels_split.all.pdf", width = 10, height = 6)
DotPlot(meninges.all.names,
        features = rev(c("C1qb","Pf4","Fscn1","Ms4a1","Iglc2","Vpreb3","Cd3g","Il7r","Gzma","Cd209a","Cox6a2","S100a9","Mki67","Col1a1","Igfbp3","Vwf","Rgs5","Enpp2","Chgb","Mpz","Plac8")),
        assay = "RNA", 
        split.by = 'orig.ident',
        cols = c("grey80", "red"), col.min = 0) + FontSize(y.title = 12, x.text = 12, y.text =12) +RotatedAxis() + xlab("") +
  ylab("Cluster Number") + scale_color_gradientn(colors = rev(rainbow(64, start = 0, end =0.7)))
dev.off()

vln.test <- VlnPlot(object = meninges.all, 
                    features = c("Siglech"), idents = 18, split.by = 'orig.ident', col = "mediumaquamarine") + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")
vln.test

vln.test <- VlnPlot(object = meninges.all, 
                    features = c("Irf7"), idents = 18, split.by = 'orig.ident', col = c("mediumaquamarine","purple")) + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")
vln.test

vln.test <- VlnPlot(object = meninges.all, 
                    features = c("Tyrobp"), idents = c(1), split.by = 'orig.ident', col = c("mediumaquamarine","purple")) + theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + xlab("")
vln.test

#now prepare objects to do differential expression analysis 
meninges.all <- readRDS("meninges.all.080320.rds")
setwd("/scratch/aco5uv/scRNAseq_meninges/LukensLab-92447355/FASTQ_Generation_2020-07-21_13_33_44Z-285512227/")
View(meninges.all)
meninges.all <- readRDS("meninges.all.080320.rds")

meninges.all.names <- readRDS("meninges.all.names.080320.rds")

View(meninges.all2)
View(meninges.all2@meta.data[["nCount_RNA"]])
View(meninges.all2@meta.data)
View(meninges.all.names)


raw.data <- as.matrix(GetAssayData(meninges.all.names,slot = "counts"))
raw.data[1:20,1:5]
saveRDS(raw.data, file = "raw.counts.data.meninges.all.080520.rds")
cluster.id <- as.matrix(Idents(meninges.all.names))
head(cluster.id)
head(t(cluster.id))
View()


#### make proportions charts ####
saveRDS(cluster.id, file = "cluster.ids.meninges.all.080520.rds")
cluster.frequency.sham.tbi <- table(Idents(meninges.all.names), meninges.all.names$orig.ident)
View(cluster.frequency.sham.tbi)
proportions.table <- (prop.table(table(Idents(meninges.all.names), meninges.all.names$orig.ident), margin = 2))
View(proportions.table)
write.table(proportions.table, file = "proportions.cells.in.cluster.011120.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(cluster.frequency.sham.tbi, file = "cells.in.cluster.080520.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
proportions.sham.tbi <- read.delim("proportions.cells.in.cluster2.080520.txt",sep = "\t")
View(proportions.sham.tbi)
proportions.data.frame <- as.data.frame(proportions.table)
View(proportions.data.frame)

bar.plot<- ggplot(proportions.data.frame, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(width = .9, stat = "identity", color = "black")
bar.plot

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI//proportion.cells.in.cluster.barchart2.pdf", width = 12, height =8)
bar.plot
dev.off()

bar.plot2<- ggplot(proportions.data.frame, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(width = .5, stat = "identity", position = position_dodge(width = .5)) + theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1), text = element_text(size=15)) + scale_fill_manual(values = c("mediumaquamarine", "purple"))
bar.plot2

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI//proportion.cells.in.cluster.groupedbarchart.pdf", width = 14, height =8)
bar.plot2
dev.off()

mac.feature.plots <- (FeaturePlot(meninges.all.names, 
                                  reduction = "umap", 
                                  features = c("C1qc", "Cd68",)))
mac.feature.plots

all.feature.plots <- (FeaturePlot(meninges.all.names, 
                                  reduction = "umap", 
                                  features = c("C1qc", "Cd68","Cd19","Cd3e","Gzmb","Col1a1","Vwf","Siglech","Plvap")))
all.feature.plots


pdf(file = "/home/aco5uv/meninges_all_outputgraphs_080320//feature.plots.signaturegenes.pdf", width = 14, height =14)
all.feature.plots
dev.off()

####Recluster Endothelial cells####
Cluster.endo.only <- subset(x = meninges.all.names, idents = c('Endothelial Cells 1','Endothelial Cells 2'))

meninges.endo <- FindVariableFeatures(Cluster.endo.only, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.endo), 10)
top10

plot1 <- VariableFeaturePlot(meninges.endo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.meninges.endo <- rownames(meninges.endo)
meninges.endo <- ScaleData(meninges.endo, features = all.genes.meninges.endo, vars.to.regress = "percent.mt")

meninges.endo <- RunPCA(meninges.endo, features = VariableFeatures(object = meninges.endo))
print(meninges.TBI[["pca"]], dims = 1:5, nfeatures = 5)

meninges.endo <- JackStraw(meninges.endo, num.replicate = 100)
meninges.endo <- ScoreJackStraw(meninges.endo, dims = 1:20)

meninges.endo <- FindNeighbors(meninges.endo, dims = 1:15)
meninges.endo <- FindClusters(meninges.endo, resolution = 0.5)
head(Idents(meninges.endo), 5)

meninges.endo <- RunUMAP(meninges.endo, dims = 1:15)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.21.UMAP.meninges.endo.split.pdf", width = 14, height =8)
DimPlot(meninges.endo, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.2021.UMAP.meninges.endo.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 3)
dev.off()

saveRDS(meninges.endo, file = "meninges.endo.011221.rds")

meninges.endo.markers <- FindAllMarkers(meninges.endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.gene <- (meninges.endo %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))
View(top.gene)

saveRDS(meninges.all.markers, file = "meninges.endo.markers.011221.rds")
write.table(meninges.endo.markers, file = "meninges_endo_markers_011221.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

endo.feature.plots <- (FeaturePlot(meninges.endo, 
                                  reduction = "umap", 
                                  features = c("Flt4")))

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.2021.UMAP.endo.featureplot.Pdpn.pdf", width = 10, height =8)
FeaturePlot(meninges.endo, reduction = "umap", features = c("Pdpn"))
dev.off()

####Recluster T cells####
meninges.all.names <- readRDS("meninges.all.names.011220.rds")
View(meninges.all.names)
Cluster.T.only <- subset(x = meninges.all.names, idents = c('Activated T Cells','CD3+ T Cells'))

meninges.T <- FindVariableFeatures(Cluster.T.only, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.T), 10)
top10

plot1 <- VariableFeaturePlot(meninges.T)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.meninges.T <- rownames(meninges.T)
meninges.T <- ScaleData(meninges.T, features = all.genes.meninges.T, vars.to.regress = "percent.mt")

meninges.T <- RunPCA(meninges.T, features = VariableFeatures(object = meninges.T))

meninges.T <- JackStraw(meninges.T, num.replicate = 100)
meninges.T <- ScoreJackStraw(meninges.T, dims = 1:20)

meninges.T <- FindNeighbors(meninges.T, dims = 1:15)
meninges.T <- FindClusters(meninges.T, resolution = 0.5)
head(Idents(meninges.T), 5)

meninges.T <- RunUMAP(meninges.T, dims = 1:15)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.21.UMAP.meninges.T.split.pdf", width = 14, height =8)
DimPlot(meninges.T, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.2021.UMAP.meninges.T.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.T, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 3)
dev.off()

saveRDS(meninges.T, file = "meninges.T.011221.rds")

meninges.T.markers <- FindAllMarkers(meninges.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(meninges.T.markers, file = "meninges.T.markers.011221.rds")
write.table(meninges.T.markers, file = "meninges_T_markers_011221.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####scCATCH T####
clu_markers <- (findmarkergenes(object = meninges.T,
                               species = 'Mouse',
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05))

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = c('Blood','Peripheral blood','Lymph node'))

View(clu_ann)


T.feature.plots <- (FeaturePlot(meninges.T, 
                                   reduction = "umap", 
                                   features = c("Flt4")))

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.12.2021.UMAP.endo.featureplot.Pdpn.pdf", width = 10, height =8)
FeaturePlot(meninges.T, reduction = "umap", features = c("Pdpn"))
dev.off()

meninges.T[["OriginalClusters"]] <- Idents(object = meninges.T)
new.cluster.ids <- c("CD8+ T Cells","Th2 Cells","NK/NKT Cells","Dying Cells","Th17 Cells")
names(new.cluster.ids) <- levels(meninges.T)
meninges.T.names <- RenameIdents(meninges.T, new.cluster.ids)

my_levels <- c("CD8+ T Cells","Th2 Cells","Th17 Cells", "NK/NKT Cells","Dying Cells")
meninges.T.names@active.ident <- factor(x=meninges.T.names@active.ident, levels = my_levels)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.T.cluster_colors_withlabels.pdf", width = 10, height =8)
DimPlot(meninges.T.names, reduction = "umap", pt.size = 5)
dev.off()

saveRDS(meninges.T.names, file = "meninges.T.names.011821.rds")

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.T.split.withnames.pdf", width = 16, height =8)
DimPlot(meninges.T.names, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham", "TBI"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2020.UMAP.meninges.T.names.WITHREPEL.pdf", width = 14, height =8)
DimPlot(meninges.T.names, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)
dev.off()

T.feature.plots <- (FeaturePlot(meninges.T.names, 
                                  reduction = "umap", 
                                  pt.size = 3,
                                  features = c("Ncr1")))
T.feature.plots

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.31.2021.Ncr1.feature.pdf", width = 8, height =8)
T.feature.plots
dev.off()

####Recluster B Cells####
meninges.all.names <- readRDS("meninges.all.names.011220.rds")
View(meninges.all.names)
View(meninges.all.names$active.ident)
Cluster.B.only <- subset(x = meninges.all.names, idents = c('B Cells 1','B Cells 2','Immature/Differentiating B Cells'))
meninges.B <- FindVariableFeatures(Cluster.B.only, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.B), 10)
top10

plot1 <- VariableFeaturePlot(meninges.B)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.meninges.B <- rownames(meninges.B)
meninges.B <- ScaleData(meninges.B, features = all.genes.meninges.B, vars.to.regress = "percent.mt")

meninges.B <- RunPCA(meninges.B, features = VariableFeatures(object = meninges.B))

meninges.B <- JackStraw(meninges.B, num.replicate = 100)
meninges.B <- ScoreJackStraw(meninges.B, dims = 1:20)

meninges.B <- FindNeighbors(meninges.B, dims = 1:15)
meninges.B <- FindClusters(meninges.B, resolution = 0.5)
head(Idents(meninges.B), 5)

meninges.B <- RunUMAP(meninges.B, dims = 1:15)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.21.UMAP.meninges.B.split.pdf", width = 14, height =8)
DimPlot(meninges.B, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.B.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.B, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 3)
dev.off()

saveRDS(meninges.B, file = "meninges.B.011821.rds")

meninges.B.markers <- FindAllMarkers(meninges.B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(meninges.B.markers, file = "meninges.B.markers.011821.rds")
write.table(meninges.B.markers, file = "meninges_B_markers_011821.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####scCATCH B####
clu_markers_B <- (findmarkergenes(object = meninges.B,
                                species = 'Mouse',
                                cluster = 'All',
                                match_CellMatch = FALSE,
                                cancer = NULL,
                                tissue = NULL,
                                cell_min_pct = 0.25,
                                logfc = 0.25,
                                pvalue = 0.05))

clu_ann_B <- scCATCH(object = clu_markers_B$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = c('Blood','Peripheral blood','Lymph node'))

View(clu_ann_B)

meninges.B[["OriginalClusters"]] <- Idents(object = meninges.B)
new.cluster.ids <- c("Mature B Cells","Activated B Cells destined to become Ab-secreting Cells","Immature B Cells","Proliferating Cells","Dying Cells")
names(new.cluster.ids) <- levels(meninges.B)
meninges.B.names <- RenameIdents(meninges.B, new.cluster.ids)

my_levels <- c("Mature B Cells","Activated B Cells destined to become Ab-secreting Cells","Immature B Cells","Proliferating Cells","Dying Cells")
meninges.B.names@active.ident <- factor(x=meninges.B.names@active.ident, levels = my_levels)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.B.cluster_colors_withlabels.pdf", width = 13, height =8)
DimPlot(meninges.B.names, reduction = "umap", pt.size = 5)
dev.off()

saveRDS(meninges.B.names, file = "meninges.B.names.011821.rds")

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.B.split.withnames.pdf", width = 16, height =8)
DimPlot(meninges.B.names, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham", "TBI"))

B.feature.plots <- (FeaturePlot(meninges.B.names, 
                                reduction = "umap", 
                                pt.size = 3,
                                features = c("Cd37")))
B.feature.plots

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.31.2021.Cd37.feature.pdf", width = 8, height =8)
B.feature.plots
dev.off()

####Recluster Macs####
meninges.all.names <- readRDS("meninges.all.names.011220.rds")
View(meninges.all.names)
Cluster.macs.only <- subset(x = meninges.all.names, idents = c('Activated Macrophages 1','Activated Macrophages 2','Macrophages 3'))

meninges.macs <- FindVariableFeatures(Cluster.macs.only, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.macs), 10)
top10

plot1 <- VariableFeaturePlot(meninges.macs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.meninges.macs <- rownames(meninges.macs)
meninges.macs <- ScaleData(meninges.macs, features = all.genes.meninges.macs, vars.to.regress = "percent.mt")

meninges.macs <- RunPCA(meninges.macs, features = VariableFeatures(object = meninges.macs))

meninges.macs <- JackStraw(meninges.macs, num.replicate = 100)
meninges.macs <- ScoreJackStraw(meninges.macs, dims = 1:20)

meninges.macs <- FindNeighbors(meninges.macs, dims = 1:15)
meninges.macs <- FindClusters(meninges.macs, resolution = 0.5)
head(Idents(meninges.macs), 5)

meninges.macs <- RunUMAP(meninges.macs, dims = 1:15)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.21.UMAP.meninges.macs.split.pdf", width = 14, height =8)
DimPlot(meninges.macs, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.18.2021.UMAP.meninges.macs.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.macs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 3)
dev.off()

saveRDS(meninges.macs, file = "meninges.macs.011221.rds")

meninges.macs.markers <- FindAllMarkers(meninges.macs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(meninges.macs.markers, file = "meninges.macs.markers.011821.rds")
write.table(meninges.macs.markers, file = "meninges_macs_markers_011821.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####scCATCH Macs####
clu_markers_macs <- (findmarkergenes(object = meninges.macs,
                                  species = 'Mouse',
                                  cluster = 'All',
                                  match_CellMatch = FALSE,
                                  cancer = NULL,
                                  tissue = NULL,
                                  cell_min_pct = 0.25,
                                  logfc = 0.25,
                                  pvalue = 0.05))

clu_ann_mac <- scCATCH(object = clu_markers_macs$clu_markers,
                     species = 'Mouse',
                     cancer = NULL,
                     tissue = c('Blood','Peripheral blood','Lymph node'))

View(clu_ann_mac)
meninges.macs[["OriginalClusters"]] <- Idents(object = meninges.macs)
new.cluster.ids <- c("Ferritin Expressing","Anti-Inflammatory","Resolution Phase","Inflammatory 1","Inflammatory 2", "Dying Cells")
names(new.cluster.ids) <- levels(meninges.macs)
meninges.macs.names <- RenameIdents(meninges.macs, new.cluster.ids)

my_levels <- c("Ferritin Expressing","Anti-Inflammatory","Resolution Phase","Inflammatory 1","Inflammatory 2", "Dying Cells")
meninges.macs.names@active.ident <- factor(x=meninges.macs.names@active.ident, levels = my_levels)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.19.2021.UMAP.meninges.macs.cluster_colors_withlabels.pdf", width = 10, height =8)
DimPlot(meninges.macs.names, reduction = "umap", pt.size = 3)
dev.off()

saveRDS(meninges.macs.names, file = "meninges.macs.names.011921.rds")

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.19.2021.UMAP.meninges.macs.split.withnames.pdf", width = 16, height =8)
DimPlot(meninges.macs.names, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', pt.size = 3, cols = rev(c("mediumaquamarine","purple")), order = c("sham", "TBI"))
dev.off()

mac.feature.plots <- (FeaturePlot(meninges.macs.names, 
                                reduction = "umap", 
                                pt.size = 2,
                                features = c("Lgals3")))
mac.feature.plots

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.31.2021.Lgals3.feature.pdf", width = 8, height =8)
mac.feature.plots
dev.off()

####frequency plots macs####

cluster.frequency.macs <- table(Idents(meninges.macs.names), meninges.macs.names$orig.ident)
View(cluster.frequency.macs)
proportions.table.macs <- (prop.table(table(Idents(meninges.macs.names), meninges.macs.names$orig.ident), margin = 2))
View(proportions.table.macs)
write.table(proportions.table.macs, file = "proportions.macs.in.cluster.011921.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(cluster.frequency.macs, file = "macs.in.cluster.011921.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
proportions.macs <- read.delim("proportions.macs.in.cluster.011921.txt",sep = "\t")
View(proportions.macs)
proportions.data.frame <- as.data.frame(proportions.macs)
proportions.data.frame2 <- as.data.frame(proportions.table.macs)
View(proportions.data.frame)

bar.plot<- ggplot(proportions.data.frame2, aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(width = .9, stat = "identity", color = "black")
bar.plot

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI//proportion.macs.in.cluster.barchart2.pdf", width = 8, height =6)
bar.plot
dev.off()

bar.plot2<- ggplot(proportions.data.frame2, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(width = .5, stat = "identity", position = position_dodge(width = .5)) + theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1), text = element_text(size=25)) + scale_fill_manual(values = c("mediumaquamarine", "purple"))
bar.plot2

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI//proportion.macs.in.cluster.groupedbarchart.pdf", width = 13, height =8)
bar.plot2
dev.off()




####Differential expression analysis####
markers.to.plot <- c("C1qb","Pf4","Fscn1","Ms4a1","Iglc2","Vpreb3","Cd3g","Il7r","Gzma","Cd209a","Cox6a2","S100a9","Mki67","Col1a1","Igfbp3","Vwf","Rgs5","Enpp2","Chgb","Mpz","Plac8")
Dotplot.identity <- DotPlot(meninges.all.names, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "orig.ident") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//01.16.21.dotplotclusteridentity_spplitbyrx.pdf", width = 10, height =8)
Dotplot.identity
dev.off()

View(meninges.all.names)
meninges.all.names$celltype.rx <- paste(Idents(meninges.all.names), meninges.all.names$orig.ident, sep = "_")
meninges.all.names$celltype <- Idents(meninges.all.names)
Idents(meninges.all.names) <- "celltype.rx"
Fibroblast.response <- FindMarkers(meninges.all.names, ident.1 = "Fibroblasts_TBI.meninges", ident.2 = "Fibroblasts_sham.meninges", verbose = FALSE)
head(Fibroblast.response, n = 15)
write.xlsx(Fibroblast.response, 'Fibroblast_diffexp_011621.xlsx', row.names=TRUE)

Actmacs.response <- FindMarkers(meninges.all.names, ident.1 = "Activated Macrophages 1_TBI.meninges", ident.2 = "Activated Macrophages 1_sham.meninges", verbose = FALSE)
head(Actmacs.response, n = 15)

####Recluster Fibroblasts####
meninges.all.names <- readRDS("meninges.all.names.011220.rds")
View(meninges.all.names)
Cluster.fibs.only <- subset(x = meninges.all.names, idents = c('Fibroblasts'))

meninges.fibs <- FindVariableFeatures(Cluster.fibs.only, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(meninges.fibs), 10)
top10

plot1 <- VariableFeaturePlot(meninges.fibs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.meninges.fibss <- rownames(meninges.fibs)
meninges.fibs <- ScaleData(meninges.fibs, features = all.genes.meninges.fibs, vars.to.regress = "percent.mt")

meninges.fibs <- RunPCA(meninges.fibs, features = VariableFeatures(object = meninges.fibs))

meninges.fibs <- JackStraw(meninges.fibs, num.replicate = 100)
meninges.fibs <- ScoreJackStraw(meninges.fibs, dims = 1:20)

meninges.fibs <- FindNeighbors(meninges.fibs, dims = 1:15)
meninges.fibs <- FindClusters(meninges.fibs, resolution = 0.5)
head(Idents(meninges.fibs), 5)

meninges.fibs <- RunUMAP(meninges.fibs, dims = 1:15)

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//02.01.21.UMAP.meninges.fibs.split.pdf", width = 14, height =8)
DimPlot(meninges.fibs, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = rev(c("mediumaquamarine","purple")), order = c("sham.meninges", "TBI.meninges"))
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/allmeninges_outputgraphs_122120//02.01.2021.UMAP.meninges.fibs.cluster_colors_WITHREPEL.pdf", width = 10, height =8)
DimPlot(meninges.fibs, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 3)
dev.off()

saveRDS(meninges.fibs, file = "meninges.fibs.020121.rds")

meninges.fibs.markers <- FindAllMarkers(meninges.fibs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(meninges.fibs.markers, file = "meninges.fibs.markers.020121.rds")
write.table(meninges.fibs.markers, file = "meninges_fibs_markers_020121.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####bulkseq and scRNAseq comparisons####
library(devtools)
library(ComplexHeatmap)
library(limma)
library(nichenetr)
library(DESeq2)
library(IRanges)
library(Seurat)
library(tidyverse)
library(gdata)

setwd("/Users/ashleybolte/Documents/Graduate School/Lukens Lab/Data/scRNAseq/scRNAseq_bulkRNAseq_connections/")
compare_bulk_sc <- read.csv("scRNA_bulkRNA_upanddown_082622.csv")
View(compare_bulk_sc_dataframe)

#finding intersections between macs, fibs, T and B cells and the upregulated genes unique to the YT vs. AT bulk commparison
intersect.YTATupunique.scRNAmacs.up<- intersect(compare_bulk_sc$Unique_YT_AT_UP,compare_bulk_sc$scRNA_macs_up)
intersect.YTATupunique.scRNAmacs.up<- as.data.frame(intersect.YTATupunique.scRNAmacs.up)
View(intersect.YTATupunique.scRNAmacs.up)

intersect.YTATupunique.scRNAfibs.up<- intersect(compare_bulk_sc$Unique_YT_AT_UP,compare_bulk_sc$scRNA_fibroblasts_up)
intersect.YTATupunique.scRNAfibs.up<- as.data.frame(intersect.YTATupunique.scRNAfibs.up)
View(intersect.YTATupunique.scRNAfibs.up)

intersect.YTATupunique.scRNABcells.up<- intersect(compare_bulk_sc$Unique_YT_AT_UP,compare_bulk_sc$scRNA_Bcell_up)
intersect.YTATupunique.scRNABcells.up<- as.data.frame(intersect.YTATupunique.scRNABcells.up)
View(intersect.YTATupunique.scRNABcells.up)

intersect.YTATupunique.scRNATcells.up<- intersect(compare_bulk_sc$Unique_YT_AT_UP,compare_bulk_sc$scRNA_Tcell_up)
intersect.YTATupunique.scRNATcells.up<- as.data.frame(intersect.YTATupunique.scRNATcells.up)
View(intersect.YTATupunique.scRNATcells.up)

total.YTATupunique <- cbindX(intersect.YTATupunique.scRNAmacs.up, intersect.YTATupunique.scRNAfibs.up, intersect.YTATupunique.scRNABcells.up,intersect.YTATupunique.scRNATcells.up)
View(total.YTATupunique)
total.YTATupunique <- total.YTATupunique[-c(1), ]
View(total.YTATupunique)
write_csv(total.YTATupunique, file = 'total_YTATupunique_082722.csv')

#finding intersections between macs, fibs, T and B cells and the downregulated genes unique to the YT vs. AT bulk commparison
intersect.YTATdownunique.scRNAmacs.down<- intersect(compare_bulk_sc$Unique_YT_AT_DOWN,compare_bulk_sc$scRNA_macs_down)
intersect.YTATdownunique.scRNAmacs.down<- as.data.frame(intersect.YTATdownunique.scRNAmacs.down)
View(intersect.YTATdownunique.scRNAmacs.down)

intersect.YTATdownunique.scRNAfibs.down<- intersect(compare_bulk_sc$Unique_YT_AT_DOWN,compare_bulk_sc$scRNA_fibroblasts_down)
intersect.YTATdownunique.scRNAfibs.down<- as.data.frame(intersect.YTATdownunique.scRNAfibs.down)
View(intersect.YTATdownunique.scRNAfibs.down)

intersect.YTATdownunique.scRNABcells.down<- intersect(compare_bulk_sc$Unique_YT_AT_DOWN,compare_bulk_sc$scRNA_Bcell_down)
intersect.YTATdownunique.scRNABcells.down<- as.data.frame(intersect.YTATdownunique.scRNABcells.down)
View(intersect.YTATdownunique.scRNABcells.down)

intersect.YTATdownunique.scRNATcells.down<- intersect(compare_bulk_sc$Unique_YT_AT_DOWN,compare_bulk_sc$scRNA_Tcell_down)
intersect.YTATdownunique.scRNATcells.down<- as.data.frame(intersect.YTATdownunique.scRNATcells.down)
View(intersect.YTATdownunique.scRNATcells.down)

total.YTATdownunique <- cbindX(intersect.YTATdownunique.scRNAmacs.down, intersect.YTATdownunique.scRNAfibs.down, intersect.YTATdownunique.scRNABcells.down,intersect.YTATdownunique.scRNATcells.down)
View(total.YTATdownunique)

total.YTATdownunique <- total.YTATdownunique[-c(1), ]
View(total.YTATdownunique)

write_csv(total.YTATdownunique, file = 'total_YTATdownunique_082722.csv')

#finding intersections between macs, fibs, T and B cells and the upregulated genes in the YS vs. AS bulk commparison
intersect.YSASup.scRNAmacs.up<- intersect(compare_bulk_sc$YS_AS_UP,compare_bulk_sc$scRNA_macs_up)
intersect.YSASup.scRNAmacs.up<- as.data.frame(intersect.YSASup.scRNAmacs.up)
View(intersect.YSASup.scRNAmacs.up)

intersect.YSASup.scRNAfibs.up<- intersect(compare_bulk_sc$YS_AS_UP,compare_bulk_sc$scRNA_fibroblasts_up)
intersect.YSASup.scRNAfibs.up<- as.data.frame(intersect.YSASup.scRNAfibs.up)
View(intersect.YSASup.scRNAfibs.up)

intersect.YSASup.scRNABcells.up<- intersect(compare_bulk_sc$YS_AS_UP,compare_bulk_sc$scRNA_Bcell_up)
intersect.YSASup.scRNABcells.up<- as.data.frame(intersect.YSASup.scRNABcells.up)
View(intersect.YSASup.scRNABcells.up)

intersect.YSASup.scRNATcells.up<- intersect(compare_bulk_sc$YS_AS_UP,compare_bulk_sc$scRNA_Tcell_up)
intersect.YSASup.scRNATcells.up<- as.data.frame(intersect.YSASup.scRNATcells.up)
View(intersect.YSASup.scRNATcells.up)

total.YSASup <- cbindX(intersect.YSASup.scRNAmacs.up, intersect.YSASup.scRNAfibs.up, intersect.YSASup.scRNABcells.up,intersect.YSASup.scRNATcells.up)
View(total.YSASup)
write_csv(total.YSASup, file = 'total_YSASup_082722.csv')

#finding intersections between macs, fibs, T and B cells and the downregulated genes in the YS vs. AS bulk commparison
intersect.YSASdown.scRNAmacs.down<- intersect(compare_bulk_sc$YS_AS_DOWN,compare_bulk_sc$scRNA_macs_down)
intersect.YSASdown.scRNAmacs.down<- as.data.frame(intersect.YSASdown.scRNAmacs.down)
View(intersect.YSASdown.scRNAmacs.down)

intersect.YSASdown.scRNAfibs.down<- intersect(compare_bulk_sc$YS_AS_DOWN,compare_bulk_sc$scRNA_fibroblasts_down)
intersect.YSASdown.scRNAfibs.down<- as.data.frame(intersect.YSASdown.scRNAfibs.down)
View(intersect.YSASdown.scRNAfibs.down)

intersect.YSASdown.scRNABcells.down<- intersect(compare_bulk_sc$YS_AS_DOWN,compare_bulk_sc$scRNA_Bcell_down)
intersect.YSASdown.scRNABcells.down<- as.data.frame(intersect.YSASdown.scRNABcells.down)
View(intersect.YSASdown.scRNABcells.down)

intersect.YSASdown.scRNATcells.down<- intersect(compare_bulk_sc$YS_AS_DOWN,compare_bulk_sc$scRNA_Tcell_down)
intersect.YSASdown.scRNATcells.down<- as.data.frame(intersect.YSASdown.scRNATcells.down)
View(intersect.YSASdown.scRNATcells.down)

total.YSASdown <- cbindX(intersect.YSASdown.scRNAmacs.down, intersect.YSASdown.scRNAfibs.down, intersect.YSASdown.scRNABcells.down,intersect.YSASdown.scRNATcells.down)
View(total.YSASdown)

write_csv(total.YSASdown, file = 'total_YSASdown_082722.csv')

####NICHENET####
library(nichenetr)
library(Seurat)
library(tidyverse)
library(viridis)

setwd("/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/")

#load in my object 
seuratobj <- readRDS("meninges.all.names.011220.rds")
View(seuratobj)
seuratobj@meta.data %>% head()
seuratobj@active.ident %>% table()
seuratobj$CellType <- Idents(seuratobj) #adding active.ident, which includes my popuation labels, into the metadata 
View(seuratobj)
seuratobj@meta.data$CellType %>% table()
DimPlot(seuratobj, reduction = "umap")
seuratobj@meta.data$orig.ident %>% table()
DimPlot(seuratobj, reduction = "umap", group.by = "orig.ident")

####long version-USE - Receiver T, Sender immune####
#load in the nichenet matrixes 
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig)

head(weighted_networks$gr)

#convert to mouse stuff 
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

#run nichenet 
#step by step version 
receiver = c("Activated T Cells","CD3+ T Cells")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3","B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_Tcellreceiver_immunesender_102422.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3","B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Dotplot_Tcell_immunecell_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_Tcellreceiver_immunesender_102422.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_regulatorypotential_Tcellreceiver_immunecellsender_Nichenet.pdf", width = 6, height =6)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetwork_Tcellreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetworkSTRICT_Tcellreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_LFCofligandsinsendercells_Tcellreceiver_immunecellsender_Nichenet.pdf", width = 4, height =6)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_pligandpearsonWITHSCALE_Tcellreceiver_immunecellsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()


#### Wrapper version NOT GOOD ####
#load in my object 
seuratobj <- readRDS("meninges.all.names.011220.rds")
View(seuratobj)
seuratobj@meta.data %>% head()
seuratobj@active.ident %>% table()
seuratobj$CellType <- Idents(seuratobj) #adding active.ident, which includes my popuation labels, into the metadata 
View(seuratobj)
seuratobj@meta.data$CellType %>% table()
DimPlot(seuratobj, reduction = "umap")
seuratobj@meta.data$orig.ident %>% table()
DimPlot(seuratobj, reduction = "umap", group.by = "orig.ident")

#load in nichenet stuff
options(timeout=100)
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig)

head(weighted_networks$gr)

#perform analysis
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratobj, 
  receiver = c("Activated T Cells","CD3+ T Cells"),
  condition_colname = "orig.ident", condition_oi = "TBI.meninges", condition_reference = "sham.meninges", 
  sender = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3","B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
View(nichenet_output)

ligand_activities = nichenet_output$ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

DotPlot(seuratobj %>% subset(idents = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3","B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

nichenet_output$ligand_differential_expression_heatmap

# Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network


#### long version- receiver B cells, sender immune ####
#run nichenet 
#step by step version 
receiver = c("B Cells 1", "B Cells 2", "Immature/Differentiating B Cells")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3","NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_Bcellreceiver_immunesender_102422.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3", "NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Dotplot_Bcellreceiver_immunecellsender_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_Bcellreceiver_immunesender_102422.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_regulatorypotential_Bcellreceiver_immunecellsender_Nichenet.pdf", width = 6, height =6)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetwork_Bcellreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetworkSTRICT_Bcellreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_LFCofligandsinsendercells_Bcellreceiver_immunecellsender_Nichenet.pdf", width = 4, height =6)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_pligandpearsonWITHSCALE_Bcellreceiver_immunecellsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()

#### long version- receiver macrophages, sender immune ####
#run nichenet 
#step by step version 
receiver = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "Fibroblasts", "Endothelial Cells 1", "Endothelial Cells 2")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_macreceiver_immunesender_102422.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "Fibroblasts", "Endothelial Cells 1", "Endothelial Cells 2")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Dotplot_macreceiver_immunecellsender_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_macreceiver_immunesender_102422.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_regulatorypotential_macreceiver_immunecellsender_Nichenet.pdf", width = 6, height =6)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetwork_macreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetworkSTRICT_macreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_LFCofligandsinsendercells_macreceiver_immunecellsender_Nichenet.pdf", width = 4, height =6)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_pligandpearsonWITHSCALE_macreceiver_immunecellsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()

#### long version- receiver fibroblasts, sender immune ####
#run nichenet 
#step by step version 
#step by step version 
receiver = c("Fibroblasts")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3", "Endothelial Cells 1", "Endothelial Cells 2")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_fibreceiver_immunesender_102422.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3", "Endothelial Cells 1", "Endothelial Cells 2")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Dotplot_fibreceiver_immunecellsender_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_fibreceiver_immunesender_102422.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_regulatorypotential_fibreceiver_immunecellsender_Nichenet.pdf", width = 6, height =6)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetwork_fibreceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_ligandreceptornetworkSTRICT_fibeceiver_immunecellsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_LFCofligandsinsendercells_fibreceiver_immunecellsender_Nichenet.pdf", width = 4, height =6)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_pligandpearsonWITHSCALE_fibreceiver_immunecellsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.24.22.Heatmap_pligandpearsonWITHoutSCALE_fibreceiver_immunecellsender_Nichenet.pdf", width = 2, height = 4)
p_ligand_pearson
dev.off()

#### long version- all cells, sender macrophages ####
#run nichenet 
#step by step version 
receiver = c("NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "Endothelial Cells 1", "Endothelial Cells 2", "Fibroblasts")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_allcellreceiver_macsender_102522.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "Endothelial Cells 1", "Endothelial Cells 2", "Fibroblasts")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Dotplot_allcellreceiver_macsender_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_allcellreceiver_macsender_102522.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090), option = "turbo")
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_regulatorypotential_allcellreceiver_macsender_Nichenet.pdf", width = 14, height =4)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_ligandreceptornetwork_allcellreceiver_macsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_ligandreceptornetworkSTRICT_allcellreceiver_macsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_LFCofligandsinsendercells_allcellreceiver_macsender_Nichenet.pdf", width = 2, height =8)
p_ligand_lfc
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_LFCofligandsinsendercells_WITSCALE_allcellreceiver_macsender_Nichenet.pdf", width = 6, height =8)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_pligandpearsonWITHSCALE_allcellreceiver_macsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_pligandpearsonWITHoutSCALE_allcellreceiver_macsender_Nichenet.pdf", width = 2, height = 12)
p_ligand_pearson
dev.off()

#### lv- all cells receiver, sender fibroblasts ####
#run nichenet 
#step by step version 
receiver = c("NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3", "Endothelial Cells 1", "Endothelial Cells 2")
expressed_genes_receiver = get_expressed_genes(receiver, seuratobj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Fibroblasts")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratobj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#define a geneset of interest 
seurat_obj_receiver= subset(seuratobj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["orig.ident"]])

condition_oi = "TBI.meninges"
condition_reference = "sham.meninges" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.1 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
View(ligand_activities)
write_csv(ligand_activities, file = 'ligandactivities_allcellreceiver_fibsender_102522.csv')

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)
Nichenet_bestligand_dotplot <- DotPlot(seuratobj %>% subset(idents = c("B Cells 1", "B Cells 2", "Immature/Differentiating B Cells", "NK Cells", "Dendritic Cells", "Activated T Cells","CD3+ T Cells", "Endothelial Cells 1", "Endothelial Cells 2", "Activated Macrophages 1","Activated Macrophages 2", "Macrophages 3")), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Dotplot_allcellreceiver_fibsender_Nichenet_bestligands.pdf", width = 10, height =4)
Nichenet_bestligand_dotplot
dev.off()

#active target gene
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.10)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

View(active_ligand_target_links)
df_active_ligand_target_links <- as.data.frame(active_ligand_target_links)
write_csv(df_active_ligand_target_links, file = 'activeligandtargetlinks_allcellreceiver_fibsender_102522.csv')

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "blue",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_viridis(breaks = c(0,0.0045,0.0090), option = "turbo")
p_ligand_target_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_regulatorypotential_allcellreceiver_fibsender_Nichenet.pdf", width = 10, height =4)
p_ligand_target_network
dev.off()

#receptors of top ranked ligands 
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
View(best_upstream_ligands)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_ligandreceptornetwork_allcellreceiver_fibsender_Nichenet.pdf", width = 8, height =6)
p_ligand_receptor_network
dev.off()

#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "DarkBlue", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_ligandreceptornetworkSTRICT_allcellreceiver_fibsender_Nichenet.pdf", width = 8, height =4)
p_ligand_receptor_network_strict
dev.off()

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratobj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratobj, condition_colname = "orig.ident", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "darkblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_LFCofligandsinsendercells_allcellreceiver_fibsender_Nichenet.pdf", width = 2, height =10)
p_ligand_lfc
dev.off()

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "DarkBlue",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_pligandpearsonWITHSCALE_allcellreceiver_fibsender_Nichenet.pdf", width = 6, height =8)
p_ligand_pearson
dev.off()

pdf(file = "/Volumes/Ashley_Backup_Lab_2/scRNAseq_Meninges_1wkpTBI/Nichenet_output_102422///10.25.22.Heatmap_pligandpearsonWITHoutSCALE_allcellreceiver_fibsender_Nichenet.pdf", width = 2, height = 12)
p_ligand_pearson
dev.off()