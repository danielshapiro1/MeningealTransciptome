#Calculate some quality control measures 
setwd("/Volumes/Ashley_Backup_Lab_2/bulkRNA_meninges_122120//")
library(lattice)
library(DESeq2) 
library(pheatmap) 
library(GSA) 
library(seq2pathway)
library(ggplot2)
library(dplyr)
library(tidyverse)

qc <- read.table('qc_metrics.txt',sep='\t',header=TRUE)
View(qc)

rownames(qc) = c('Aged Sham 1','Aged Sham 2','Aged Sham 3','Aged TBI 1', 'Aged TBI 2','Aged TBI 3', 'Young Sham 1','Young Sham 2','Young Sham 3', 'Young TBI 1','Young TBI 2','Young TBI 3')

qc = qc[,-1]

qc$prop.junk = 1 - (qc$without.junk / qc$all.counts)
qc$prop.dups = 1 - (qc$unique.counts / qc$without.junk)

load('all.counts.Rdata') #this object isn't created until later 
qc$prop.in.features = colSums(all.counts[1:55421,]) / colSums(all.counts)

qc$colors = c(rep('dark grey',3),rep('orange',3),rep('light blue',3),rep('blue',3))


pdf('total.reads.pdf', width=10, height=4)
print(barchart(all.counts ~ rownames(qc), data = qc,
               main = "Total Reads",
               xlab = "Sample",
               ylab = "Read Count",
               col = qc$colors,
               ylim = c(0,100000000),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.junk.pdf', width=10, height=4)
print(barchart(prop.junk ~ rownames(qc), data = qc,
               main = "Proportion Junk",
               xlab = "Sample",
               ylab = "Proportion Junk",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.duplicates.pdf', width=10, height=4)
print(barchart(prop.dups ~ rownames(qc), data = qc,
               main = "Proportion PCR Duplicates",
               xlab = "Sample",
               ylab = "Proportion PCR Duplicates",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.in.features.pdf', width=10, height=4)
print(barchart(prop.in.features ~ rownames(qc), data = qc,
               main = "Proportion In Features",
               xlab = "Sample",
               ylab = "Proportion In Features",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

#we create a function to define our genes as not changed, or up or downregulated, this is with padj
categorize.deseq.df <- function(df, thresh = 0.1, log2fold = 0.0, treat = 'Auxin') {
  df.activated = data.frame(matrix(nrow = 0, ncol = 0))
  df.repressed = data.frame(matrix(nrow = 0, ncol = 0))
  if (nrow(df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold,]) != 0) {
    df.activated = df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold,]
    df.activated$response = paste(treat, 'Activated')
  }
  if (nrow(df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold,]) != 0) {
    df.repressed = df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold,]
    df.repressed$response = paste(treat, 'Repressed')
  }
  df.unchanged = df[df$padj > 0.5 & !is.na(df$padj) & abs(df$log2FoldChange) < 0.25,]
  df.unchanged$response = paste(treat, 'Unchanged')
  df.dregs = df[!(df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold) &
                  !(df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold) &
                  !(df$padj > 0.5 & !is.na(df$padj) &
                      abs(df$log2FoldChange) < 0.25), ]
  df.dregs$response = paste(treat, 'All Other Genes')
  df.effects.lattice =
    rbind(df.activated,
          df.unchanged,
          df.repressed,
          df.dregs)
  df.effects.lattice$response = factor(df.effects.lattice$response)
  df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'Unchanged'))
  df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'All Other Genes'))
  return(df.effects.lattice)
}

#this is with pval 
categorize.deseq.df2 <- function(df, thresh = 0.05, log2fold = 0.0, treat = 'Auxin') {
  df.activated = data.frame(matrix(nrow = 0, ncol = 0))
  df.repressed = data.frame(matrix(nrow = 0, ncol = 0))
  if (nrow(df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold,]) != 0) {
    df.activated = df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold,]
    df.activated$response = paste(treat, 'Activated')
  }
  if (nrow(df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold,]) != 0) {
    df.repressed = df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold,]
    df.repressed$response = paste(treat, 'Repressed')
  }
  df.unchanged = df[df$pval > 0.5 & !is.na(df$pval) & abs(df$log2FoldChange) < 0.25,]
  df.unchanged$response = paste(treat, 'Unchanged')
  df.dregs = df[!(df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold) &
                  !(df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold) &
                  !(df$pval > 0.5 & !is.na(df$pval) &
                      abs(df$log2FoldChange) < 0.25), ]
  df.dregs$response = paste(treat, 'All Other Genes')
  df.effects.lattice2 =
    rbind(df.activated,
          df.unchanged,
          df.repressed,
          df.dregs)
  df.effects.lattice2$response = factor(df.effects.lattice2$response)
  df.effects.lattice2$response = relevel(df.effects.lattice2$response, ref = paste(treat, 'Unchanged'))
  df.effects.lattice2$response = relevel(df.effects.lattice2$response, ref = paste(treat, 'All Other Genes'))
  return(df.effects.lattice2)
}

#PCA plot function
plotPCAlattice <- function(df, file = 'PCA_lattice.pdf') {
  perVar = round(100 * attr(df, "percentVar"))
  df = data.frame(cbind(df, sapply(strsplit(as.character(df$name), '_rep'), '[', 1)))
  colnames(df) = c(colnames(df)[1:(ncol(df)-1)], 'unique_condition')
  print(df)
  color.x = substring(rainbow(length(unique(df$unique_condition))), 1,7)
  df$color = NA
  df$alpha.x = NA
  df$alpha.y = NA
  df$colpal = NA
  for (i in 1:length(unique(df$unique_condition))) {
    df[df$unique_condition == unique(df$unique_condition)[[i]],]$color = color.x[i]
    reps_col<- df[df$unique_condition == unique(df$unique_condition)[[i]],]
    replicates.x = nrow(reps_col)
    alx <- rev(seq(0.2, 1, length.out = replicates.x))
    for(rep in 1:replicates.x) {
      na <- reps_col[rep, ]$name
      df[df$name == na, ]$alpha.x = alx[rep]
      aly = as.hexmode(round(alx * 255))
      df[df$name == na, ]$alpha.y = aly[rep]
      cp = paste0(color.x[i], aly)
      df[df$name == na, ]$colpal = cp[rep]
    }
  }
  colpal = df$colpal
  df$name = gsub('_', ' ', df$name)
  df$name <- factor(df$name, levels=df$name, order=TRUE)
  pdf(file, width=6, height=6, useDingbats=FALSE)
  print(xyplot(PC2 ~ PC1, groups = name, data=df,
               xlab = paste('PC1: ', perVar[1], '% variance', sep = ''),
               ylab = paste('PC2: ', perVar[2], '% variance', sep = ''),
               par.settings = list(superpose.symbol = list(pch = c(20), col=colpal)),
               pch = 20, cex = 1.7,
               aspect = 1,
               auto.key = TRUE,
               col = colpal))
  dev.off()
}

#now we are loading the data in and creating the dataframe
x <- read.table("Aged_Sham_rep1.gene.counts.txt",header=FALSE)

all.counts <- data.frame(row.names = x$V1)

samples <- c("Aged_Sham_rep1", "Aged_Sham_rep2","Aged_Sham_rep3","Aged_TBI_rep1", "Aged_TBI_rep2","Aged_TBI_rep3","Young_Sham_rep1","Young_Sham_rep2","Young_Sham_rep3","Young_TBI_rep1","Young_TBI_rep2","Young_TBI_rep3")

for (i in samples){
  j<-read.table(print(paste0(i,".gene.counts.txt")))
  all.counts <- cbind(all.counts,data.frame(i = j[2]))
}

colnames(all.counts) <- samples

save(all.counts,file='all.counts.Rdata')
View(all.counts)

tail(all.counts)
merged.counts <- all.counts[-(55386:55390),]
View(merged.counts)
tail(merged.counts)

ensembl.all = read.table('gencode.vM24.annotation.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');

save(ensembl.code, file = "ensembl.code.Rdata")

rm(ensembl.all)
rm(ensembl.gene.names)
rm(ensembl.id.names)

merged.counts = merge(merged.counts, ensembl.code, by="row.names", all.x=F)
rownames(merged.counts) <- make.names(merged.counts$gene, unique=TRUE)
merged.counts <- merged.counts[,-c(1,14,15)]

rm(ensembl.code)
View(merged.counts)

save(merged.counts, file="merged.counts.Rdata")

unt=3
trt=9
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)
pca.dat = plotPCA(rld_HH, intgroup="condition", returnData=TRUE)
percentVar = round(100 * attr(pca.dat, "percentVar"))
plotPCAlattice(pca.dat, file = 'PCA_RNAseq.pdf')

####Differential expression analysis and creating MA plots#####

####YoungTBI vs. AgedTBI####
merged.counts.small = merged.counts[,c(10:12,4:6)]
View(merged.counts.small)

unt = 3
trt = 3

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.youngTBIvagedTBI <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

young.tbi.v.aged.tbi.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'Aged')

pdf("young.tbi.v.aged.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(young.tbi.v.aged.tbi.lattice$log2FoldChange ~
               log(young.tbi.v.aged.tbi.lattice$baseMean, base=10),
             groups=young.tbi.v.aged.tbi.lattice$response,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+Age log"[2]~"Expression fold change"),
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'Young_TBI v. Aged_TBI',
             par.settings=list(par.xlab.text=list(cex=1.1,font=2),
                               par.ylab.text=list(cex=1.1,font=2), strip.background=list(col="grey85"))))
dev.off()

save(young.tbi.v.aged.tbi.lattice, file='young.tbi.v.aged.tbi.lattice.Rdata')

####Aged_Sham v Aged_TBI####
merged.counts.small = merged.counts[,1:6]
View(merged.counts.small)

unt = 3
trt = 3

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

aged.sham.v.aged.tbi.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'TBI')

pdf("aged.sham.v.aged.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(aged.sham.v.aged.tbi.lattice$log2FoldChange ~
               log(aged.sham.v.aged.tbi.lattice$baseMean, base=10),
             groups=aged.sham.v.aged.tbi.lattice$response,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+TBI log"[2]~"Expression fold change"),
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'Aged Sham v. Aged TBI',
             par.settings=list(par.xlab.text=list(cex=1.1,font=2),
                               par.ylab.text=list(cex=1.1,font=2),
                               strip.background=list(col="grey85"))))
dev.off()

save(aged.sham.v.aged.tbi.lattice, file='aged.sham.v.aged.tbi.lattice.Rdata')

####Young_Sham vs Young_TBI####
View(merged.counts)
merged.counts.small = merged.counts[,c(7:9,10:12)]
View(merged.counts.small)

unt = 3
trt = 3

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)


mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)


sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)


young.sham.v.young.tbi.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'TBI')

pdf("young.sham.v.young.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(young.sham.v.young.tbi.lattice$log2FoldChange ~
               log(young.sham.v.young.tbi.lattice$baseMean, base=10),
             groups=young.sham.v.young.tbi.lattice$response,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+TBI log"[2]~"Expression fold change"),
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'Young Sham v. Young TBI',
             par.settings=list(par.xlab.text=list(cex=1.1,font=2),
                               par.ylab.text=list(cex=1.1,font=2),
                               strip.background=list(col="grey85"))))
dev.off()

save(young.sham.v.young.tbi.lattice, file='young.sham.v.young.tbi.lattice.Rdata')

####Young_Sham v Aged_Sham####
View(merged.counts)
merged.counts.small = merged.counts[,c(7:9,1:3)]
View(merged.counts.small)

unt = 3
trt = 3

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

young.sham.v.aged.sham.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'Age')

pdf("young.sham.v.aged.sham.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(young.sham.v.aged.sham.lattice$log2FoldChange ~
               log(young.sham.v.aged.sham.lattice$baseMean, base=10),
             groups=young.sham.v.aged.sham.lattice$response,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+Age log"[2]~"Expression fold change"),
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'Young Sham v. Aged Sham',
             par.settings=list(par.xlab.text=list(cex=1.1,font=2),
                               par.ylab.text=list(cex=1.1,font=2),
                               strip.background=list(col="grey85"))))
dev.off()

save(young.sham.v.aged.sham.lattice, file='young.sham.v.aged.sham.lattice.Rdata')

#### # Differentially Expressed Genes ####
#show the number of differentially expressed genes
print(table(young.tbi.v.aged.tbi.lattice$response))
print(table(aged.sham.v.aged.tbi.lattice$response))
print(table(young.sham.v.young.tbi.lattice$response))
print(table(young.sham.v.aged.sham.lattice$response))

####table of genes####
#make a table of the increased and decreased genes for GO terms
lattice = young.tbi.v.aged.tbi.lattice

write(rownames(lattice[lattice$response == 'Aged Activated',]),file='0.1.padj.inc.txt',sep='\t')
write(rownames(lattice[lattice$response == 'Aged Repressed',]),file='0.1.padj.dec.txt',sep='\t')

lattice = aged.sham.v.aged.tbi.lattice

write(rownames(lattice[lattice$response == 'TBI Activated',]),file='aged.sham.v.aged.tbi.0.1.padj.inc.txt',sep='\t')
write(rownames(lattice[lattice$response == 'TBI Repressed',]),file='aged.sham.v.aged.tbi.0.1.padj.dec.txt',sep='\t')

lattice = young.sham.v.young.tbi.lattice

write(rownames(lattice[lattice$response == 'TBI Activated',]),file='young.sham.v.young.tbi.0.1.padj.inc.txt',sep='\t')
write(rownames(lattice[lattice$response == 'TBI Repressed',]),file='young.sham.v.young.tbi.0.1.padj.dec.txt',sep='\t')

lattice = young.sham.v.aged.sham.lattice

write(rownames(lattice[lattice$response == 'Aged Activated',]),file='young.sham.v.aged.sham.0.1.padj.inc.txt',sep='\t')
write(rownames(lattice[lattice$response == 'Aged Repressed',]),file='young.sham.v.aged.sham.0.1.padj.dec.txt',sep='\t')

#### young v aged tbi heatmap####
#we will make heatmaps now
unt=3
trt=3

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts[,c(10:12,4:6)], DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = young.tbi.v.aged.tbi.lattice
lattice = na.omit(lattice)

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:20]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:20]

a = assay(rld_HH)[c(inc[1:20],dec[1:20]),]
a = a - rowMeans(a)
View(a)

pdf(file='young.v.aged.tbi.top20.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=T))
dev.off()

#### young v aged sham heatmap####
#we will make heatmaps now
unt=3
trt=3
View(merged.counts)
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts[,c(7:9,1:3)], DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = young.sham.v.aged.sham.lattice
lattice = na.omit(lattice)

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:20]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:20]

a = assay(rld_HH)[c(inc[1:20],dec[1:20]),]
a = a - rowMeans(a)
View(a)

pdf(file='young.v.aged.sham.top20.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=T))
dev.off()

#### aged sham v aged tbi heatmap####
unt=3
trt=3
View(merged.counts)
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts[,c(1:3,4:6)], DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = aged.sham.v.aged.tbi.lattice
lattice = na.omit(lattice)

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:20]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:20]

a = assay(rld_HH)[c(inc[1:20],dec[1:20]),]
a = a - rowMeans(a)
View(a)

pdf(file='aged.sham.v.aged.tbi.top20.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=T))
dev.off()

####volcano plots####
volcano.plot <- function(lattice, thresh = 0.1, treat='Aged', title="", highlights = c()) {
  name = strsplit(deparse(substitute(lattice)),'.lattice')[[1]][1]
  
  lattice$color = 'black'
  lattice[lattice$response == paste0(treat,' Activated'),]$color = 'red'
  lattice[lattice$response == paste0(treat,' Repressed'),]$color = 'blue'
  
  if (length(highlights) > 0) {
    name = paste0(name,'.highlights')
    highlights = paste0('^',highlights,'$')
    row.nums = c()
    for(i in highlights) {
      row.nums = append(row.nums, grep(i,rownames(lattice))) }
    lattice[row.nums,]$color = 'green'
  }
  
  pdf(file = paste0(name,'.volcano.plot.pdf'))
  
  plot(lattice$log2FoldChange, -log10(lattice$padj), xlim = c(-5,6), col = lattice$color, cex = 0.8, pch = 20,
       main = title, xlab = 'log2FoldChange', ylab = '-log10(padj)')
  
  abline(h=-log(thresh,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
  dev.off()
}

View(young.tbi.v.aged.tbi.lattice)
volcano.plot(young.tbi.v.aged.tbi.lattice, treat='Aged', title='Young TBI v. Aged TBI')
volcano.plot(aged.sham.v.aged.tbi.lattice, treat='TBI', title='Aged Sham v. Aged TBI')
volcano.plot(young.sham.v.young.tbi.lattice, treat ='TBI', title='Young Sham v. Young TBI')
volcano.plot(young.sham.v.aged.sham.lattice, treat='Age', title='Young Sham v. Aged Sham')

View(young.sham.v.aged.sham.lattice)

lattice = young.sham.v.aged.sham.lattice
lattice
active.aged.sham <- as.data.frame((lattice[lattice$response == 'Age Activated',c(2,5,6)]))
View(active.aged)
library(openxlsx)
write.xlsx(active.aged, 'young.sham.v.aged.sham.aged.activated.xlsx', row.names=TRUE)

write.xlsx(merged.counts, 'merged.counts.122920.xlsx', row.names=TRUE)

####make plot from GO terms####
GOterms.young.sham.v.aged.sham <- read.csv("go_biological_process_youngsham_v_agedsham.csv")
View(GOterms.young.sham.v.aged.sham)

#filter dataset 
filtered.GO<- (GOterms.young.sham.v.aged.sham %>% 
                     arrange(adjusted_p_value))
filtered.GO<- (GOterms.young.sham.v.aged.sham %>% top_n(-25, adjusted_p_value))
View(filtered.GO)

#graph data 
young.v.age.sham <- (ggplot(filtered.GO, 
                    aes(x = negative_log10_of_adjusted_p_value, y =fct_reorder(term_name, -adjusted_p_value))) + 
               geom_point(aes(size = term_size, color = intersection_size)) +
               theme_bw(base_size = 14) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
               scale_colour_gradient(limits=c(90, 500), low="blue", high="green") +
               ylab(NULL) +
               ggtitle("Young Sham v. Aged Sham Enriched GO Terms"))
young.v.age.sham
ggsave(
  path = "/Volumes/Ashley_Backup_Lab_2/bulkRNA_meninges_122120/",
  filename = "youngsham.v.agedsham.GO.pdf",
  plot = young.v.age.sham, 
  device = "pdf", 
  dpi=150, 
  width=185, 
  height=185, 
  units = "mm")

####find intersections between Young Sham v Aged Sham and Young TBI v. Aged TBI ####

lattice = young.sham.v.aged.sham.lattice
active.aged.sham <- as.data.frame((lattice[lattice$response == 'Age Activated',c(2,5,6)]))

lattice = young.tbi.v.aged.tbi.lattice
active.aged.tbi <- as.data.frame((lattice[lattice$response == 'Aged Activated',c(2,5,6)]))
View(active.aged.tbi)

intersections.activated <- intersect((rownames(active.aged.tbi)),(rownames(active.aged.sham)))
unique.activated <- unique((rownames(active.aged.tbi)),(intersections.activated))
View(intersections.activated)
View(active.aged.tbi)
View(active.aged.sham)
View(unique.activated)

View(young.sham.v.aged.sham.lattice)

lattice = young.sham.v.aged.sham.lattice
depressed.aged.sham <- as.data.frame((lattice[lattice$response == 'Age Repressed',c(2,5,6)]))
View(depressed.aged.sham)

lattice = young.tbi.v.aged.tbi.lattice
depressed.aged.tbi <- as.data.frame((lattice[lattice$response == 'Aged Repressed',c(2,5,6)]))
View(depressed.aged.tbi)

intersections.repressed <- intersect((rownames(depressed.aged.sham)),(rownames(depressed.aged.tbi)))
unique.repressed <- unique((rownames(active.aged.tbi)),(intersections.activated))
View(intersections.repressed)

setwd("/Volumes/Ashley_Backup_Lab_2/bulkRNA_meninges_122120//")
View(intersections.activated)
View(intersections.repressed)
typeof(intersections.repressed)
all.intersections <- (intersections.activated %>% unite(intersections.repressed, sep = "", remove = FALSE)) 
write.xlsx(intersections.activated, 'intersectionsactivated.xlsx', row.names=TRUE)
write.xlsx(intersections.repressed, 'intersectionsrepressed.xlsx', row.names=TRUE)
intersection.all <- (read_csv('intersectionsall.csv', col_names = FALSE))
View(intersection.all)

lattice = young.tbi.v.aged.tbi.lattice
allchanged.tbi <- as.data.frame((lattice[lattice$response %in% c('Aged Repressed','Aged Activated'),c(2,5,6,7)]))
View(allchanged.tbi)
library(data.table)
allchanged.tbi.row <- rownames_to_column(allchanged.tbi, var = "genes") %>% as_tibble()
View(allchanged.tbi.row)
allchanged.tbi.onlynames <- (allchanged.tbi.row[,1])
View(allchanged.tbi.onlynames)
allchanged.tbi.onlynames <- as_vector((allchanged.tbi.onlynames))
intersection.all <- as_vector(intersection.all)
differences.tbi <- setdiff(allchanged.tbi.onlynames,intersection.all)
View(differences.tbi)
write_csv(differences.tbi.dataframe, file = 'differences.tbi.123120')
View(allchanged.tbi.row)
differences.tbi.dataframe <- as_data_frame(differences.tbi)
View(differences.tbi.dataframe)
datawithoutrows = allchanged.tbi[ !(rownames(allchanged.tbi) %in% intersection.all), ]
View(datawithoutrows)
library(openxlsx)
write.xlsx(datawithoutrows, 'YoungTBIvsAgedTBI.minus.YoungShamvsAgedSham.xlsx', row.names=TRUE)

####plot out stuff unique to TBI And Aging####
####repressed####
GOterms.repressed.agedtbi <- read.csv("repressed_uniquetoAgedTBI_gomolecularfunction_123120.csv")
View(GOterms.repressed.agedtbi)

#filter dataset 
filtered.GO.repressed<- (GOterms.repressed.agedtbi %>% 
                 arrange(adjusted_p_value))
filtered.GO.repressed<- (GOterms.repressed.agedtbi %>% top_n(-25, adjusted_p_value))
View(filtered.GO.repressed)

#graph data 
aged.tbi.repressed <- (ggplot(filtered.GO.repressed, 
                            aes(x = negative_log10_of_adjusted_p_value, y =fct_reorder(term_name, -adjusted_p_value))) + 
                       geom_point(aes(size = term_size, color = intersection_size, )) +
                       theme_bw(base_size = 14) +
                       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                       scale_colour_gradient(limits=c(10, 902), low="blue", high="green") +
                       ylab(NULL) +
                       ggtitle("Unique Repressed Genes in Aged TBI"))
aged.tbi.repressed
ggsave(
  path = "/Volumes/Ashley_Backup_Lab_2/bulkRNA_meninges_122120/",
  filename = "Unique_repressed_genes_Aged_TBI.pdf",
  plot = aged.tbi.repressed, 
  device = "pdf", 
  dpi=150, 
  width=185, 
  height=185, 
  units = "mm")

####activated####
GOterms.activated.agedtbi <- read.csv("Go_biologicalProcess_uniquetoAgedTBI_010121.csv")
View(GOterms.activated.agedtbi)

#filter dataset 
filtered.GO.activated<- (GOterms.activated.agedtbi %>% 
                           arrange(adjusted_p_value))
filtered.GO.activated<- (GOterms.activated.agedtbi %>% top_n(-10, adjusted_p_value))
View(filtered.GO.activated)

#graph data 
aged.tbi.activated <- (ggplot(filtered.GO.activated, 
                              aes(x = fct_reorder(term_name, -adjusted_p_value), y =negative_log10_of_adjusted_p_value)) + 
                         geom_bar(aes(fill = intersection_size), stat = 'identity') +
                         theme_bw(base_size = 14) +
                         theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 1)) + 
                         scale_fill_viridis(option = "viridis") +
                         ylab(NULL) +
                         ggtitle("Unique Activated Genes in Aged TBI"))
aged.tbi.activated
ggsave(
  path = "/Volumes/Ashley_Backup_Lab_2/bulkRNA_meninges_122120/",
  filename = "Unique_activated_genes_Aged_TBI_bar.pdf",
  plot = aged.tbi.activated, 
  device = "pdf", 
  dpi=150, 
  width=200, 
  height=175, 
  units = "mm")
