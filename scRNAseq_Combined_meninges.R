library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)

setwd("/Users/danielshapiro/Desktop/Research/Code/")
sham.data <- Read10X(data.dir = "ShamData")
sham <- CreateSeuratObject(counts = sham.data, project = "sham.meninges", min.cells = 3, min.features = 200)
TBI.data <- Read10X(data.dir = "TBIData")
TBI <- CreateSeuratObject(counts = TBI.data, project = "TBI.meninges", min.cells = 3, min.features = 200)

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


pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.QCgraphs.meninges.all.pdf")
(qc.featurerna.meninges.all | qc.nCountrna.meninges.all | qc.percentMt.meninges.all | qc.percentHemo.meninges.all) / (qc.count.featurerna.scatter.meninges.all | qc.count.percentmt.scatter.meninges.all | qc.count.percenthemo.scatter.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.featurena.meninges.all.pdf")
(qc.featurerna.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.nCountrna.meninges.all.pdf")
(qc.nCountrna.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.percentMt.meninges.all.pdf")
(qc.percentMt.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.percentHemo.meninges.all.pdf")
(qc.percentHemo.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.count.featurerna.scatter.meninges.all.pdf")
(qc.count.featurerna.scatter.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.count.percentmt.scatter.meninges.all.pdf")
(qc.count.percentmt.scatter.meninges.all)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.qc.count.percenthemo.scatter.meninges.all.pdf")
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

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.variablefeatures.meninges.all.pdf")
(plot2)
dev.off()

all.genes.meninges.all <- rownames(meninges.all)
meninges.all <- ScaleData(meninges.all, features = all.genes.meninges.all, vars.to.regress = "percent.mt")

meninges.all <- RunPCA(meninges.all, features = VariableFeatures(object = meninges.all))
print(meninges.all[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.Vizdimreduction.pca.meninges.all.pdf")
VizDimLoadings(meninges.all, dims = 1:2, reduction = "pca")
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs//08.03.2020.dimplot.pca.meninges.all.pdf")
DimPlot(meninges.all, reduction = "pca")
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.Dimheatmap.10dim.meninges.all.pdf")
DimHeatmap(meninges.all, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.Dimheatmap.1dim.meninges.all.pdf")
DimHeatmap(meninges.all, dims = 1, cells = 500, balanced = TRUE)
dev.off()

meninges.all <- JackStraw(meninges.all, num.replicate = 100)
meninges.all <- ScoreJackStraw(meninges.all, dims = 1:20)

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.jackstrawdiminsionality.all.pdf")
JackStrawPlot(meninges.all, dims = 1:20)
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.01.2020.elbowplotdiminsionality.all.pdf")
ElbowPlot(meninges.all)
dev.off()

meninges.all <- FindNeighbors(meninges.all, dims = 1:15)
meninges.all <- FindClusters(meninges.all, resolution = 0.5)
head(Idents(meninges.all), 5)

meninges.all <- RunUMAP(meninges.all, dims = 1:15)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', cols = c("cornflowerblue", "darkseagreen2"), order = c("TBI", "Sham"))

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.pdf", width = 10, height =8)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', cols = c("darkturquoise", "orchid1"), order = c("TBI", "Sham"))
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.split.pdf", width = 14, height =8)
DimPlot(meninges.all, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = c("darkturquoise", "orchid1"), order = c("TBI", "Sham"))
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.cluster_colors.pdf", width = 10, height =8)
DimPlot(meninges.all, reduction = "umap", label = TRUE)
dev.off()

View(meninges.all)
saveRDS(meninges.all, file = "meninges.all.080320.rds")

meninges.all.markers <- FindAllMarkers(meninges.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.gene <- (meninges.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))
View(top.gene)

#---

saveRDS(meninges.all.markers, file = "meninges.all.markers.080320.rds")
write.table(meninges.all.markers, file = "meninges_all_markers_080120.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

meninges.all[["OriginalClusters"]] <- Idents(object = meninges.all)
new.cluster.ids <- c("Fibroblasts","B Cells","Activated Macrophages","Vascular Endothelial Cells","Activated Macrophages 2","CD3+ T Cells","Dendritic Cells","Activated T cells","Endothelial Cells","Immature/Differentiating B cells","Ciliated Ependymal Cells/Choroid Plexus (Pia mater)","NK Cells","Pineal Gland Cells","Pericytes","Schwann/Neural Crest Cells","Macrophages 3","Plasmacytoid Dendritic Cells/IFN response","Macrophages 4","Neutrophils")
names(new.cluster.ids) <- levels(meninges.all)
meninges.all.names <- RenameIdents(meninges.all, new.cluster.ids)

my_levels <- c("Fibroblasts","B Cells","Activated Macrophages","Vascular Endothelial Cells","Activated Macrophages 2","CD3+ T Cells","Dendritic Cells","Activated T cells","Endothelial Cells","Immature/Differentiating B cells","Ciliated Ependymal Cells/Choroid Plexus (Pia mater)","NK Cells","Pineal Gland Cells","Pericytes","Schwann/Neural Crest Cells","Macrophages 3","Plasmacytoid Dendritic Cells/IFN response","Macrophages 4","Neutrophils")
meninges.all.names@active.ident <- factor(x=meninges.all.names@active.ident, levels = my_levels)

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.cluster_colors_withlabels.pdf", width = 14, height =8)
DimPlot(meninges.all.names, reduction = "umap")
dev.off()

saveRDS(meninges.all.names, file = "meninges.all.names.080320.rds")

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.split.withnames.pdf", width = 16, height =8)
DimPlot(meninges.all.names, reduction = "umap", group.by ='orig.ident', split.by ='orig.ident', cols = c("darkturquoise", "orchid1"), order = c("TBI", "Sham"))
dev.off()

meninges.all.names[["OriginalClusters"]] <- Idents(object = meninges.all.names)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
names(new.cluster.ids) <- levels(meninges.all.names)
meninges.all.numbers <- RenameIdents(meninges.all.names, new.cluster.ids)

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.UMAP.meninges.all.cluster_colors_withnumbers.pdf", width = 10, height =8)
DimPlot(meninges.all.numbers, reduction = "umap", label = TRUE)
dev.off()

#find the markers of each cluster 
meninges.all.markers.names <- FindAllMarkers(meninges.all.names, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top.gene.names.all <- (meninges.all.markers.names %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))
View(top.gene.names.all)
write.table(top.gene.names.all, file = "top_gene_names_bycluster_all_080320.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)


#make dot plot

print(top.gene.names.all[,7])

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.dotplot_topgenes_withlabels.all.pdf", width = 10, height = 6)
DotPlot(meninges.all.names,
        features = rev(c("Col1a1","Cd79a","C1qb","Igfbp3","Pf4","Cd3e","Cd209a","Il7r","Vwf","Vpreb3","Enpp2","Gzma","Chgb","Rgs5","Mpz","Plac8","Cox6a2","Fscn1","S100a9")),
       assay = "RNA", 
        cols = c("grey80", "red"), col.min = 0) + FontSize(y.title = 12, x.text = 12, y.text =12) +RotatedAxis() + xlab("") +
  ylab("Cluster Number") + scale_color_gradientn(colors = rev(rainbow(64, start = 0, end =0.7)))
dev.off()

pdf(file = "/home/ds3vz/OutputGraphs/08.03.2020.dotplot_topgenes_withlabels_split.all.pdf", width = 10, height = 6)
DotPlot(meninges.all.names,
        features = rev(c("Col1a1","Cd79a","C1qb","Igfbp3","Pf4","Cd3e","Cd209a","Il7r","Vwf","Vpreb3","Enpp2","Gzma","Chgb","Rgs5","Mpz","Plac8","Cox6a2","Fscn1","S100a9")),
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
setwd("/scratch/ds3vz/LukensLab-92447355/FASTQ_Generation_2020-07-21_13_33_44Z-285512227/")
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

meninges.all.names <- readRDS("cluster.ids.meninges.all.080520.rds")
saveRDS(cluster.id, file = "cluster.ids.meninges.all.080520.rds")
cluster.frequency.sham.tbi <- table(Idents(meninges.all.names), meninges.all.names$orig.ident)
View(cluster.frequency.sham.tbi)
proportions.table <- (prop.table(table(Idents(meninges.all.names), meninges.all.names$orig.ident), margin = 2))
View(proportions.table)
write.table(proportions.table, file = "proportions.cells.in.cluster.080520.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(cluster.frequency.sham.tbi, file = "cells.in.cluster.080520.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
proportions.sham.tbi <- read.delim("proportions.cells.in.cluster2.080520.txt",sep = "\t")
View(proportions.sham.tbi)
proportions.data.frame <- as.data.frame(proportions.table)
View(proportions.data.frame)
a <- prop.test(x=c(127,784),n=c(5697,5697))

stat_pvale_manual()
bar.plot<- ggplot(proportions.data.frame, aes(x=Var2, y=Freq, fill=Var1))
bar.plot

pdf(file = "/home/ds3vz/OutputGraphs/proportion.cells.in.cluster.barchart2.pdf", width = 12, height =8)
bar.plot
dev.off()

bar.plot2<- ggplot(proportions.data.frame, aes(x=Var1, y=Freq, fill=Var2))+ stat_pvalue_manual(a,label="Pvalue")
  geom_bar(width = .5, stat = "identity", position = position_dodge(width = .5)) + theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1), text = element_text(size=15)) + scale_fill_manual(values = c("darkturquoise", "orchid1"))
bar.plot2

pdf(file = "/home/ds3vz/OutputGraphs/proportion.cells.in.cluster.groupedbarchart.pdf", width = 14, height =8)
bar.plot2
dev.off()

prop.test()



mac.feature.plots <- (FeaturePlot(meninges.all.names, 
            reduction = "umap", 
            features = c("C1qc", "Cd68")))
mac.feature.plots

all.feature.plots <- (FeaturePlot(meninges.all.names, 
                                  reduction = "umap", 
                                  features = c("C1qc", "Cd68","Cd19","Cd3e","Gzmb","Col1a1","Vwf","Siglech","Plvap")))
all.feature.plots


pdf(file = "/home/ds3vz/OutputGraphs/feature.plots.signaturegenes.pdf", width = 14, height =14)
all.feature.plots
dev.off()
