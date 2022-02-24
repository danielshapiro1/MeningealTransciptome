library(zinbwave)
library(DESeq2)
library(scran)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(writexl)
setwd("/Users/danielshapiro/Desktop/Code/")
categorize.deseq.df <- function(df, fdr = 0.05, log2fold = 0.0, treat
= 'Auxin') {

    df.effects.lattice = df
    df.effects.lattice$response = 'X'

    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange > log2fold,]$response = 'up'
    }
    if(nrow(df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]) > 0) {
        df.effects.lattice[df.effects.lattice$padj < fdr & !is.na(df.effects.lattice$padj) & df.effects.lattice$log2FoldChange < log2fold,]$response = 'down'
    }
     
    return(df.effects.lattice)
}

clusters = readRDS('cluster.ids.meninges.all.080520.rds')
counts = readRDS('raw.counts.data.meninges.all.080520.rds')

#### Activated Fibroblasts ####
act.tcells = clusters[clusters[,1] == "Activated T cells",]
small.counts.act.tcells = counts[,colnames(counts) %in% names(act.tcells)]
small.counts.act.tcells = small.counts.act.tcells[rowSums(small.counts.act.tcells) > 0,]

filter.act.tcells <- rowSums(small.counts.act.tcells >= 5) >= 5
table(filter.act.tcells)
small.counts.act.tcells <- small.counts.act.tcells[filter.act.tcells,]
dim(small.counts.act.tcells)

nSham.act.tcells = length(grep('Sham',colnames(small.counts.act.tcells)))
nTBI.act.tcells = length(grep('TBI',colnames(small.counts.act.tcells)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.act.tcells <- DataFrame(Treatment=c(rep('Sham',nSham.act.tcells),rep('TBI',nTBI.act.tcells)),
                     row.names=colnames(small.counts.act.tcells))

x.act.tcells = SummarizedExperiment(assays=list(counts=small.counts.act.tcells), colData=colData.act.tcells)

zinb.act.tcells <- zinbwave(x.act.tcells, K=0, epsilon=1000)

save(zinb.act.tcells,file='zinb.act.tcells.Rdata')

dds.act.tcells <- DESeqDataSet(zinb.act.tcells, design = ~ Treatment)
dds.act.tcells <- estimateSizeFactors(dds.act.tcells, type="poscounts")
#scr <- computeSumFactors(dds.act.tcells)
#sizeFactors(dds) <- sizeFactors(scr)

dds.act.tcells <- DESeq(dds.act.tcells, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.act.tcells)
res.act.tcells <- results(dds.act.tcells, independentFiltering=FALSE)
summary(res.act.tcells)
write.csv(as.data.frame(results(dds.act.tcells)),file="DESeq/act.tcells.csv")

act.t.cell.deseq = categorize.deseq.df(res.act.tcells,treat='TBI')
table(act.t.cell.deseq$response)


#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Activated Macrophages ####
act.macs = clusters[clusters[,1] == c("Activated Macrophages","Activated Macrophages 2"),]
small.counts.act.macs = counts[,colnames(counts) %in% names(act.macs)]
small.counts.act.macs = small.counts.act.macs[rowSums(small.counts.act.macs) > 0,]

filter.act.macs <- rowSums(small.counts.act.macs >= 5) >= 5
table(filter.act.macs)
small.counts.act.macs <- small.counts.act.macs[filter.act.macs,]
dim(small.counts.act.macs)

nSham.act.macs = length(grep('Sham',colnames(small.counts.act.macs)))
nTBI.act.macs = length(grep('TBI',colnames(small.counts.act.macs)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.act.macs <- DataFrame(Treatment=c(rep('Sham',nSham.act.macs),rep('TBI',nTBI.act.macs)),
                                row.names=colnames(small.counts.act.macs))

x.act.macs = SummarizedExperiment(assays=list(counts=small.counts.act.macs), colData=colData.act.macs)

zinb.act.macs <- zinbwave(x.act.macs, K=0, epsilon=1000)

save(zinb.act.macs,file='zinb.act.macs.Rdata')

dds.act.macs <- DESeqDataSet(zinb.act.macs, design = ~ Treatment)
dds.act.macs <- estimateSizeFactors(dds.act.macs, type="poscounts")
#scr <- computeSumFactors(dds.act.macs)
#sizeFactors(dds) <- sizeFactors(scr)

dds.act.macs <- DESeq(dds.act.macs, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.act.macs)
res.act.macs <- results(dds.act.macs, independentFiltering=FALSE)
summary(res.act.macs)
write.csv(as.data.frame(results(dds.act.macs)),file="DESeq/act.macs.csv")

act.t.cell.deseq = categorize.deseq.df(res.act.macs,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### B Cells ####
bcells = clusters[clusters[,1] == "B Cells",]
small.counts.bcells = counts[,colnames(counts) %in% names(bcells)]
small.counts.bcells = small.counts.bcells[rowSums(small.counts.bcells) > 0,]

filter.bcells <- rowSums(small.counts.bcells >= 5) >= 5
table(filter.bcells)
small.counts.bcells <- small.counts.bcells[filter.bcells,]
dim(small.counts.bcells)

nSham.bcells = length(grep('Sham',colnames(small.counts.bcells)))
nTBI.bcells = length(grep('TBI',colnames(small.counts.bcells)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.bcells <- DataFrame(Treatment=c(rep('Sham',nSham.bcells),rep('TBI',nTBI.bcells)),
                              row.names=colnames(small.counts.bcells))

x.bcells = SummarizedExperiment(assays=list(counts=small.counts.bcells), colData=colData.bcells)

zinb.bcells <- zinbwave(x.bcells, K=0, epsilon=1000)

save(zinb.bcells,file='zinb.bcells.Rdata')

dds.bcells <- DESeqDataSet(zinb.bcells, design = ~ Treatment)
dds.bcells <- estimateSizeFactors(dds.bcells, type="poscounts")
#scr <- computeSumFactors(dds.bcells)
#sizeFactors(dds) <- sizeFactors(scr)

dds.bcells <- DESeq(dds.bcells, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.bcells)
res.bcells <- results(dds.bcells, independentFiltering=FALSE)
summary(res.bcells)
write.csv(as.data.frame(results(dds.bcells)),file="DESeq/bcells.csv")

act.t.cell.deseq = categorize.deseq.df(res.bcells,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### CD3+ T Cells ####
cd3tcells = clusters[clusters[,1] == "CD3+ T Cells",]
small.counts.cd3tcells = counts[,colnames(counts) %in% names(cd3tcells)]
small.counts.cd3tcells = small.counts.cd3tcells[rowSums(small.counts.cd3tcells) > 0,]

filter.cd3tcells <- rowSums(small.counts.cd3tcells >= 5) >= 5
table(filter.cd3tcells)
small.counts.cd3tcells <- small.counts.cd3tcells[filter.cd3tcells,]
dim(small.counts.cd3tcells)

nSham.cd3tcells = length(grep('Sham',colnames(small.counts.cd3tcells)))
nTBI.cd3tcells = length(grep('TBI',colnames(small.counts.cd3tcells)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.cd3tcells <- DataFrame(Treatment=c(rep('Sham',nSham.cd3tcells),rep('TBI',nTBI.cd3tcells)),
                            row.names=colnames(small.counts.cd3tcells))

x.cd3tcells = SummarizedExperiment(assays=list(counts=small.counts.cd3tcells), colData=colData.cd3tcells)

zinb.cd3tcells <- zinbwave(x.cd3tcells, K=0, epsilon=1000)

save(zinb.cd3tcells,file='zinb.cd3tcells.Rdata')

dds.cd3tcells <- DESeqDataSet(zinb.cd3tcells, design = ~ Treatment)
dds.cd3tcells <- estimateSizeFactors(dds.cd3tcells, type="poscounts")
#scr <- computeSumFactors(dds.cd3tcells)
#sizeFactors(dds) <- sizeFactors(scr)

dds.cd3tcells <- DESeq(dds.cd3tcells, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.cd3tcells)
res.cd3tcells <- results(dds.cd3tcells, independentFiltering=FALSE)
summary(res.cd3tcells)
write.csv(as.data.frame(results(dds.cd3tcells)),file="DESeq/cd3tcells.csv")

act.t.cell.deseq = categorize.deseq.df(res.cd3tcells,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Ciliated Ependymal Cells/Choroid Plexus (Pia mater) ####
ependymal = clusters[clusters[,1] == "Ciliated Ependymal Cells/Choroid Plexus (Pia mater)",]
small.counts.ependymal = counts[,colnames(counts) %in% names(ependymal)]
small.counts.ependymal = small.counts.ependymal[rowSums(small.counts.ependymal) > 0,]

filter.ependymal <- rowSums(small.counts.ependymal >= 5) >= 5
table(filter.ependymal)
small.counts.ependymal <- small.counts.ependymal[filter.ependymal,]
dim(small.counts.ependymal)

nSham.ependymal = length(grep('Sham',colnames(small.counts.ependymal)))
nTBI.ependymal = length(grep('TBI',colnames(small.counts.ependymal)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.ependymal <- DataFrame(Treatment=c(rep('Sham',nSham.ependymal),rep('TBI',nTBI.ependymal)),
                               row.names=colnames(small.counts.ependymal))

x.ependymal = SummarizedExperiment(assays=list(counts=small.counts.ependymal), colData=colData.ependymal)

zinb.ependymal <- zinbwave(x.ependymal, K=0, epsilon=1000)

save(zinb.ependymal,file='zinb.ependymal.Rdata')

dds.ependymal <- DESeqDataSet(zinb.ependymal, design = ~ Treatment)
dds.ependymal <- estimateSizeFactors(dds.ependymal, type="poscounts")
#scr <- computeSumFactors(dds.ependymal)
#sizeFactors(dds) <- sizeFactors(scr)

dds.ependymal <- DESeq(dds.ependymal, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.ependymal)
res.ependymal <- results(dds.ependymal, independentFiltering=FALSE)
summary(res.ependymal)
write.csv(as.data.frame(results(dds.ependymal)),file="DESeq/ependymal.csv")

act.t.cell.deseq = categorize.deseq.df(res.ependymal,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Vascular Endothelial Cells ####
vasc.endo = clusters[clusters[,1] == "Vascular Endothelial Cells",]
small.counts.vasc.endo = counts[,colnames(counts) %in% names(vasc.endo)]
small.counts.vasc.endo = small.counts.vasc.endo[rowSums(small.counts.vasc.endo) > 0,]

filter.vasc.endo <- rowSums(small.counts.vasc.endo >= 5) >= 5
table(filter.vasc.endo)
small.counts.vasc.endo <- small.counts.vasc.endo[filter.vasc.endo,]
dim(small.counts.vasc.endo)

nSham.vasc.endo = length(grep('Sham',colnames(small.counts.vasc.endo)))
nTBI.vasc.endo = length(grep('TBI',colnames(small.counts.vasc.endo)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.vasc.endo <- DataFrame(Treatment=c(rep('Sham',nSham.vasc.endo),rep('TBI',nTBI.vasc.endo)),
                               row.names=colnames(small.counts.vasc.endo))

x.vasc.endo = SummarizedExperiment(assays=list(counts=small.counts.vasc.endo), colData=colData.vasc.endo)

zinb.vasc.endo <- zinbwave(x.vasc.endo, K=0, epsilon=1000)

save(zinb.vasc.endo,file='zinb.vasc.endo.Rdata')

dds.vasc.endo <- DESeqDataSet(zinb.vasc.endo, design = ~ Treatment)
dds.vasc.endo <- estimateSizeFactors(dds.vasc.endo, type="poscounts")
#scr <- computeSumFactors(dds.vasc.endo)
#sizeFactors(dds) <- sizeFactors(scr)

dds.vasc.endo <- DESeq(dds.vasc.endo, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.vasc.endo)
res.vasc.endo <- results(dds.vasc.endo, independentFiltering=FALSE)
summary(res.vasc.endo)
write.csv(as.data.frame(results(dds.vasc.endo)),file="DESeq/vasc.endo.csv")

act.t.cell.deseq = categorize.deseq.df(res.vasc.endo,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Dendritic Cells ####
dendritic = clusters[clusters[,1] == "Dendritic Cells",]
small.counts.dendritic = counts[,colnames(counts) %in% names(dendritic)]
small.counts.dendritic = small.counts.dendritic[rowSums(small.counts.dendritic) > 0,]

filter.dendritic <- rowSums(small.counts.dendritic >= 5) >= 5
table(filter.dendritic)
small.counts.dendritic <- small.counts.dendritic[filter.dendritic,]
dim(small.counts.dendritic)

nSham.dendritic = length(grep('Sham',colnames(small.counts.dendritic)))
nTBI.dendritic = length(grep('TBI',colnames(small.counts.dendritic)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.dendritic <- DataFrame(Treatment=c(rep('Sham',nSham.dendritic),rep('TBI',nTBI.dendritic)),
                               row.names=colnames(small.counts.dendritic))

x.dendritic = SummarizedExperiment(assays=list(counts=small.counts.dendritic), colData=colData.dendritic)

zinb.dendritic <- zinbwave(x.dendritic, K=0, epsilon=1000)

save(zinb.dendritic,file='zinb.dendritic.Rdata')

dds.dendritic <- DESeqDataSet(zinb.dendritic, design = ~ Treatment)
dds.dendritic <- estimateSizeFactors(dds.dendritic, type="poscounts")
#scr <- computeSumFactors(dds.dendritic)
#sizeFactors(dds) <- sizeFactors(scr)

dds.dendritic <- DESeq(dds.dendritic, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.dendritic)
res.dendritic <- results(dds.dendritic, independentFiltering=FALSE)
summary(res.dendritic)
write.csv(as.data.frame(results(dds.dendritic)),file="DESeq/dendritic.csv")

act.t.cell.deseq = categorize.deseq.df(res.dendritic,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Plasmacytoid Dendritic Cells/IFN response ####
plasmacytoid.dendritic = clusters[clusters[,1] == "Plasmacytoid Dendritic Cells/IFN response",]
small.counts.plasmacytoid.dendritic = counts[,colnames(counts) %in% names(plasmacytoid.dendritic)]
small.counts.plasmacytoid.dendritic = small.counts.plasmacytoid.dendritic[rowSums(small.counts.plasmacytoid.dendritic) > 0,]

filter.plasmacytoid.dendritic <- rowSums(small.counts.plasmacytoid.dendritic >= 5) >= 5
table(filter.plasmacytoid.dendritic)
small.counts.plasmacytoid.dendritic <- small.counts.plasmacytoid.dendritic[filter.plasmacytoid.dendritic,]
dim(small.counts.plasmacytoid.dendritic)

nSham.plasmacytoid.dendritic = length(grep('Sham',colnames(small.counts.plasmacytoid.dendritic)))
nTBI.plasmacytoid.dendritic = length(grep('TBI',colnames(small.counts.plasmacytoid.dendritic)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.plasmacytoid.dendritic <- DataFrame(Treatment=c(rep('Sham',nSham.plasmacytoid.dendritic),rep('TBI',nTBI.plasmacytoid.dendritic)),
                               row.names=colnames(small.counts.plasmacytoid.dendritic))

x.plasmacytoid.dendritic = SummarizedExperiment(assays=list(counts=small.counts.plasmacytoid.dendritic), colData=colData.plasmacytoid.dendritic)

zinb.plasmacytoid.dendritic <- zinbwave(x.plasmacytoid.dendritic, K=0, epsilon=1000)

save(zinb.plasmacytoid.dendritic,file='zinb.plasmacytoid.dendritic.Rdata')

dds.plasmacytoid.dendritic <- DESeqDataSet(zinb.plasmacytoid.dendritic, design = ~ Treatment)
dds.plasmacytoid.dendritic <- estimateSizeFactors(dds.plasmacytoid.dendritic, type="poscounts")
#scr <- computeSumFactors(dds.plasmacytoid.dendritic)
#sizeFactors(dds) <- sizeFactors(scr)

dds.plasmacytoid.dendritic <- DESeq(dds.plasmacytoid.dendritic, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.plasmacytoid.dendritic)
res.plasmacytoid.dendritic <- results(dds.plasmacytoid.dendritic, independentFiltering=FALSE)
summary(res.plasmacytoid.dendritic)
write.csv(as.data.frame(results(dds.plasmacytoid.dendritic)),file="DESeq/plasmacytoid.dendritic.csv")

act.t.cell.deseq = categorize.deseq.df(res.plasmacytoid.dendritic,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Macrophages 3 ####
macro3 = clusters[clusters[,1] == "Macrophages 3",]
small.counts.macro3 = counts[,colnames(counts) %in% names(macro3)]
small.counts.macro3 = small.counts.macro3[rowSums(small.counts.macro3) > 0,]

filter.macro3 <- rowSums(small.counts.macro3 >= 5) >= 5
table(filter.macro3)
small.counts.macro3 <- small.counts.macro3[filter.macro3,]
dim(small.counts.macro3)

nSham.macro3 = length(grep('Sham',colnames(small.counts.macro3)))
nTBI.macro3 = length(grep('TBI',colnames(small.counts.macro3)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.macro3 <- DataFrame(Treatment=c(rep('Sham',nSham.macro3),rep('TBI',nTBI.macro3)),
                                            row.names=colnames(small.counts.macro3))

x.macro3 = SummarizedExperiment(assays=list(counts=small.counts.macro3), colData=colData.macro3)

zinb.macro3 <- zinbwave(x.macro3, K=0, epsilon=1000)

save(zinb.macro3,file='zinb.macro3.Rdata')

dds.macro3 <- DESeqDataSet(zinb.macro3, design = ~ Treatment)
dds.macro3 <- estimateSizeFactors(dds.macro3, type="poscounts")
#scr <- computeSumFactors(dds.macro3)
#sizeFactors(dds) <- sizeFactors(scr)

dds.macro3 <- DESeq(dds.macro3, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.macro3)
res.macro3 <- results(dds.macro3, independentFiltering=FALSE)
summary(res.macro3)
write.csv(as.data.frame(results(dds.macro3)),file="DESeq/macro3.csv")

act.t.cell.deseq = categorize.deseq.df(res.macro3,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### NK Cells ####
NKcells = clusters[clusters[,1] == "NK Cells",]
small.counts.NKcells = counts[,colnames(counts) %in% names(NKcells)]
small.counts.NKcells = small.counts.NKcells[rowSums(small.counts.NKcells) > 0,]

filter.NKcells <- rowSums(small.counts.NKcells >= 5) >= 5
table(filter.NKcells)
small.counts.NKcells <- small.counts.NKcells[filter.NKcells,]
dim(small.counts.NKcells)

nSham.NKcells = length(grep('Sham',colnames(small.counts.NKcells)))
nTBI.NKcells = length(grep('TBI',colnames(small.counts.NKcells)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.NKcells <- DataFrame(Treatment=c(rep('Sham',nSham.NKcells),rep('TBI',nTBI.NKcells)),
                            row.names=colnames(small.counts.NKcells))

x.NKcells = SummarizedExperiment(assays=list(counts=small.counts.NKcells), colData=colData.NKcells)

zinb.NKcells <- zinbwave(x.NKcells, K=0, epsilon=1000)

save(zinb.NKcells,file='zinb.NKcells.Rdata')

dds.NKcells <- DESeqDataSet(zinb.NKcells, design = ~ Treatment)
dds.NKcells <- estimateSizeFactors(dds.NKcells, type="poscounts")
#scr <- computeSumFactors(dds.NKcells)
#sizeFactors(dds) <- sizeFactors(scr)

dds.NKcells <- DESeq(dds.NKcells, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.NKcells)
res.NKcells <- results(dds.NKcells, independentFiltering=FALSE)
summary(res.NKcells)
write.csv(as.data.frame(results(dds.NKcells)),file="DESeq/NKcells.csv")

act.t.cell.deseq = categorize.deseq.df(res.NKcells,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Schwann/Neural Crest Cells ####
schwann.NC = clusters[clusters[,1] == "Schwann/Neural Crest Cells",]
small.counts.schwann.NC = counts[,colnames(counts) %in% names(schwann.NC)]
small.counts.schwann.NC = small.counts.schwann.NC[rowSums(small.counts.schwann.NC) > 0,]

filter.schwann.NC <- rowSums(small.counts.schwann.NC >= 5) >= 5
table(filter.schwann.NC)
small.counts.schwann.NC <- small.counts.schwann.NC[filter.schwann.NC,]
dim(small.counts.schwann.NC)

nSham.schwann.NC = length(grep('Sham',colnames(small.counts.schwann.NC)))
nTBI.schwann.NC = length(grep('TBI',colnames(small.counts.schwann.NC)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.schwann.NC <- DataFrame(Treatment=c(rep('Sham',nSham.schwann.NC),rep('TBI',nTBI.schwann.NC)),
                             row.names=colnames(small.counts.schwann.NC))

x.schwann.NC = SummarizedExperiment(assays=list(counts=small.counts.schwann.NC), colData=colData.schwann.NC)

zinb.schwann.NC <- zinbwave(x.schwann.NC, K=0, epsilon=1000)

save(zinb.schwann.NC,file='zinb.schwann.NC.Rdata')

dds.schwann.NC <- DESeqDataSet(zinb.schwann.NC, design = ~ Treatment)
dds.schwann.NC <- estimateSizeFactors(dds.schwann.NC, type="poscounts")
#scr <- computeSumFactors(dds.schwann.NC)
#sizeFactors(dds) <- sizeFactors(scr)

dds.schwann.NC <- DESeq(dds.schwann.NC, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.schwann.NC)
res.schwann.NC <- results(dds.schwann.NC, independentFiltering=FALSE)
summary(res.schwann.NC)
write.csv(as.data.frame(results(dds.schwann.NC)),file="DESeq/schwann.NC.csv")

act.t.cell.deseq = categorize.deseq.df(res.schwann.NC,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Fibroblasts ####
fibroblasts = clusters[clusters[,1] == "Fibroblasts",]
small.counts.fibroblasts = counts[,colnames(counts) %in% names(fibroblasts)]
small.counts.fibroblasts = small.counts.fibroblasts[rowSums(small.counts.fibroblasts) > 0,]

filter.fibroblasts <- rowSums(small.counts.fibroblasts >= 5) >= 5
table(filter.fibroblasts)
small.counts.fibroblasts <- small.counts.fibroblasts[filter.fibroblasts,]
dim(small.counts.fibroblasts)

nSham.fibroblasts = length(grep('Sham',colnames(small.counts.fibroblasts)))
nTBI.fibroblasts = length(grep('TBI',colnames(small.counts.fibroblasts)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.fibroblasts <- DataFrame(Treatment=c(rep('Sham',nSham.fibroblasts),rep('TBI',nTBI.fibroblasts)),
                                row.names=colnames(small.counts.fibroblasts))

x.fibroblasts = SummarizedExperiment(assays=list(counts=small.counts.fibroblasts), colData=colData.fibroblasts)

zinb.fibroblasts <- zinbwave(x.fibroblasts, K=0, epsilon=1000)

save(zinb.fibroblasts,file='zinb.fibroblasts.Rdata')

dds.fibroblasts <- DESeqDataSet(zinb.fibroblasts, design = ~ Treatment)
dds.fibroblasts <- estimateSizeFactors(dds.fibroblasts, type="poscounts")
#scr <- computeSumFactors(dds.fibroblasts)
#sizeFactors(dds) <- sizeFactors(scr)

dds.fibroblasts <- DESeq(dds.fibroblasts, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.fibroblasts)
res.fibroblasts <- results(dds.fibroblasts, independentFiltering=FALSE)
summary(res.fibroblasts)
write.csv(as.data.frame(results(dds.fibroblasts)),file="DESeq/fibroblasts.csv")

act.t.cell.deseq = categorize.deseq.df(res.fibroblasts,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Pineal Gland Cells ####
pineal = clusters[clusters[,1] == "Pineal Gland Cells",]
small.counts.pineal = counts[,colnames(counts) %in% names(pineal)]
small.counts.pineal = small.counts.pineal[rowSums(small.counts.pineal) > 0,]

filter.pineal <- rowSums(small.counts.pineal >= 5) >= 5
table(filter.pineal)
small.counts.pineal <- small.counts.pineal[filter.pineal,]
dim(small.counts.pineal)

nSham.pineal = length(grep('Sham',colnames(small.counts.pineal)))
nTBI.pineal = length(grep('TBI',colnames(small.counts.pineal)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.pineal <- DataFrame(Treatment=c(rep('Sham',nSham.pineal),rep('TBI',nTBI.pineal)),
                                row.names=colnames(small.counts.pineal))

x.pineal = SummarizedExperiment(assays=list(counts=small.counts.pineal), colData=colData.pineal)

zinb.pineal <- zinbwave(x.pineal, K=0, epsilon=1000)

save(zinb.pineal,file='zinb.pineal.Rdata')

dds.pineal <- DESeqDataSet(zinb.pineal, design = ~ Treatment)
dds.pineal <- estimateSizeFactors(dds.pineal, type="poscounts")
#scr <- computeSumFactors(dds.pineal)
#sizeFactors(dds) <- sizeFactors(scr)

dds.pineal <- DESeq(dds.pineal, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.pineal)
res.pineal <- results(dds.pineal, independentFiltering=FALSE)
summary(res.pineal)
write.csv(as.data.frame(results(dds.pineal)),file="DESeq/pineal.csv")

act.t.cell.deseq = categorize.deseq.df(res.pineal,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Pericytes ####
pericytes = clusters[clusters[,1] == "Pericytes",]
small.counts.pericytes = counts[,colnames(counts) %in% names(pericytes)]
small.counts.pericytes = small.counts.pericytes[rowSums(small.counts.pericytes) > 0,]

filter.pericytes <- rowSums(small.counts.pericytes >= 5) >= 5
table(filter.pericytes)
small.counts.pericytes <- small.counts.pericytes[filter.pericytes,]
dim(small.counts.pericytes)

nSham.pericytes = length(grep('Sham',colnames(small.counts.pericytes)))
nTBI.pericytes = length(grep('TBI',colnames(small.counts.pericytes)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.pericytes <- DataFrame(Treatment=c(rep('Sham',nSham.pericytes),rep('TBI',nTBI.pericytes)),
                                row.names=colnames(small.counts.pericytes))

x.pericytes = SummarizedExperiment(assays=list(counts=small.counts.pericytes), colData=colData.pericytes)

zinb.pericytes <- zinbwave(x.pericytes, K=0, epsilon=1000)

save(zinb.pericytes,file='zinb.pericytes.Rdata')

dds.pericytes <- DESeqDataSet(zinb.pericytes, design = ~ Treatment)
dds.pericytes <- estimateSizeFactors(dds.pericytes, type="poscounts")
#scr <- computeSumFactors(dds.pericytes)
#sizeFactors(dds) <- sizeFactors(scr)

dds.pericytes <- DESeq(dds.pericytes, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.pericytes)
res.pericytes <- results(dds.pericytes, independentFiltering=FALSE)
summary(res.pericytes)
write.csv(as.data.frame(results(dds.pericytes)),file="DESeq/pericytes.csv")

act.t.cell.deseq = categorize.deseq.df(res.pericytes,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Immature/Differentiating B cells ####
im.dif.Bcells = clusters[clusters[,1] == "Immature/Differentiating B cells",]
small.counts.im.dif.Bcells = counts[,colnames(counts) %in% names(im.dif.Bcells)]
small.counts.im.dif.Bcells = small.counts.im.dif.Bcells[rowSums(small.counts.im.dif.Bcells) > 0,]

filter.im.dif.Bcells <- rowSums(small.counts.im.dif.Bcells >= 5) >= 5
table(filter.im.dif.Bcells)
small.counts.im.dif.Bcells <- small.counts.im.dif.Bcells[filter.im.dif.Bcells,]
dim(small.counts.im.dif.Bcells)

nSham.im.dif.Bcells = length(grep('Sham',colnames(small.counts.im.dif.Bcells)))
nTBI.im.dif.Bcells = length(grep('TBI',colnames(small.counts.im.dif.Bcells)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.im.dif.Bcells <- DataFrame(Treatment=c(rep('Sham',nSham.im.dif.Bcells),rep('TBI',nTBI.im.dif.Bcells)),
                                row.names=colnames(small.counts.im.dif.Bcells))

x.im.dif.Bcells = SummarizedExperiment(assays=list(counts=small.counts.im.dif.Bcells), colData=colData.im.dif.Bcells)

zinb.im.dif.Bcells <- zinbwave(x.im.dif.Bcells, K=0, epsilon=1000)

save(zinb.im.dif.Bcells,file='zinb.im.dif.Bcells.Rdata')

dds.im.dif.Bcells <- DESeqDataSet(zinb.im.dif.Bcells, design = ~ Treatment)
dds.im.dif.Bcells <- estimateSizeFactors(dds.im.dif.Bcells, type="poscounts")
#scr <- computeSumFactors(dds.im.dif.Bcells)
#sizeFactors(dds) <- sizeFactors(scr)

dds.im.dif.Bcells <- DESeq(dds.im.dif.Bcells, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.im.dif.Bcells)
res.im.dif.Bcells <- results(dds.im.dif.Bcells, independentFiltering=FALSE)
summary(res.im.dif.Bcells)
write.csv(as.data.frame(results(dds.im.dif.Bcells)),file="DESeq/im.dif.Bcells.csv")

act.t.cell.deseq = categorize.deseq.df(res.im.dif.Bcells,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Endothelial Cells ####
endo = clusters[clusters[,1] == "Endothelial Cells",]
small.counts.endo = counts[,colnames(counts) %in% names(endo)]
small.counts.endo = small.counts.endo[rowSums(small.counts.endo) > 0,]

filter.endo <- rowSums(small.counts.endo >= 5) >= 5
table(filter.endo)
small.counts.endo <- small.counts.endo[filter.endo,]
dim(small.counts.endo)

nSham.endo = length(grep('Sham',colnames(small.counts.endo)))
nTBI.endo = length(grep('TBI',colnames(small.counts.endo)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.endo <- DataFrame(Treatment=c(rep('Sham',nSham.endo),rep('TBI',nTBI.endo)),
                                row.names=colnames(small.counts.endo))

x.endo = SummarizedExperiment(assays=list(counts=small.counts.endo), colData=colData.endo)

zinb.endo <- zinbwave(x.endo, K=0, epsilon=1000)

save(zinb.endo,file='zinb.endo.Rdata')

dds.endo <- DESeqDataSet(zinb.endo, design = ~ Treatment)
dds.endo <- estimateSizeFactors(dds.endo, type="poscounts")
#scr <- computeSumFactors(dds.endo)
#sizeFactors(dds) <- sizeFactors(scr)

dds.endo <- DESeq(dds.endo, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.endo)
res.endo <- results(dds.endo, independentFiltering=FALSE)
summary(res.endo)
write.csv(as.data.frame(results(dds.endo)),file="DESeq/endo.csv")

act.t.cell.deseq = categorize.deseq.df(res.endo,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### Macrophages 4 ####
macro4 = clusters[clusters[,1] == "Macrophages 4",]
small.counts.macro4 = counts[,colnames(counts) %in% names(macro4)]
small.counts.macro4 = small.counts.macro4[rowSums(small.counts.macro4) > 0,]

filter.macro4 <- rowSums(small.counts.macro4 >= 5) >= 5
table(filter.macro4)
small.counts.macro4 <- small.counts.macro4[filter.macro4,]
dim(small.counts.macro4)

nSham.macro4 = length(grep('Sham',colnames(small.counts.macro4)))
nTBI.macro4 = length(grep('TBI',colnames(small.counts.macro4)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.macro4 <- DataFrame(Treatment=c(rep('Sham',nSham.macro4),rep('TBI',nTBI.macro4)),
                                row.names=colnames(small.counts.macro4))

x.macro4 = SummarizedExperiment(assays=list(counts=small.counts.macro4), colData=colData.macro4)

zinb.macro4 <- zinbwave(x.macro4, K=0, epsilon=1000)

save(zinb.macro4,file='zinb.macro4.Rdata')

dds.macro4 <- DESeqDataSet(zinb.macro4, design = ~ Treatment)
dds.macro4 <- estimateSizeFactors(dds.macro4, type="poscounts")
#scr <- computeSumFactors(dds.macro4)
#sizeFactors(dds) <- sizeFactors(scr)

dds.macro4 <- DESeq(dds.macro4, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.macro4)
res.macro4 <- results(dds.macro4, independentFiltering=FALSE)
summary(res.macro4)
write.csv(as.data.frame(results(dds.macro4)),file="DESeq/macro4.csv")

act.t.cell.deseq = categorize.deseq.df(res.macro4,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#### Neutrophils ####
neutrophils = clusters[clusters[,1] == "Neutrophils",]
small.counts.neutrophils = counts[,colnames(counts) %in% names(neutrophils)]
small.counts.neutrophils = small.counts.neutrophils[rowSums(small.counts.neutrophils) > 0,]

filter.neutrophils <- rowSums(small.counts.neutrophils >= 5) >= 5
table(filter.neutrophils)
small.counts.neutrophils <- small.counts.neutrophils[filter.neutrophils,]
dim(small.counts.neutrophils)

nSham.neutrophils = length(grep('Sham',colnames(small.counts.neutrophils)))
nTBI.neutrophils = length(grep('TBI',colnames(small.counts.neutrophils)))

#nGenes = nrow(counts)
#coverage = colSums(counts) / nGenes

colData.neutrophils <- DataFrame(Treatment=c(rep('Sham',nSham.neutrophils),rep('TBI',nTBI.neutrophils)),
                                row.names=colnames(small.counts.neutrophils))

x.neutrophils = SummarizedExperiment(assays=list(counts=small.counts.neutrophils), colData=colData.neutrophils)

zinb.neutrophils <- zinbwave(x.neutrophils, K=0, epsilon=1000)

save(zinb.neutrophils,file='zinb.neutrophils.Rdata')

dds.neutrophils <- DESeqDataSet(zinb.neutrophils, design = ~ Treatment)
dds.neutrophils <- estimateSizeFactors(dds.neutrophils, type="poscounts")
#scr <- computeSumFactors(dds.neutrophils)
#sizeFactors(dds) <- sizeFactors(scr)

dds.neutrophils <- DESeq(dds.neutrophils, test="LRT", reduced=~1, useT=TRUE, minmu=1e-6, minRep=Inf)
results(dds.neutrophils)
res.neutrophils <- results(dds.neutrophils, independentFiltering=FALSE)
summary(res.neutrophils)
write.csv(as.data.frame(results(dds.neutrophils)),file="DESeq/neutrophils.csv")

act.t.cell.deseq = categorize.deseq.df(res.neutrophils,treat='TBI')
table(act.t.cell.deseq$response)

#Counts
par(mfrow=c(2,3))
plotCounts(dds, gene="Rpl7", intgroup="Treatment")
plotCounts(dds, gene="Ptpn18", intgroup="Treatment")
plotCounts(dds, gene="Cox5b", intgroup="Treatment")
plotCounts(dds, gene="Rpl31", intgroup="Treatment")
plotCounts(dds, gene="Il1rl1", intgroup="Treatment")
plotCounts(dds, gene="mt-Co3", intgroup="Treatment")

#Volcano Plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


