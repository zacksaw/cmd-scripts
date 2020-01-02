##Differential Analysis using featureCounts
# Differential gene expression analysis with DESeq2 #

### LIBRARIES
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("genefilter")
library("biomaRt")
library("IHW")

### SET WORKING DIRECTORY
# note: this directory should be populated with the output.txt counts file generated in the prelim step
#setwd("PATH TO WORKING DIRECTORY")
setwd("/Users/zachariasawaged/Documents/MS/bioinformatics/Assignments/Assignment_2_3/Assignment3")

### Import count table and details on experimental design
# NB: Make sure column names are exactly as in counts file (same name, same order)
CountTable <- read.table("output.txt", header=TRUE, row.names=1) 
#remove extra columns from feature counts output - this may or may not have already been done
CountTable$Chr <- NULL
CountTable$Start <- NULL
CountTable$End <- NULL
CountTable$Strand <- NULL
CountTable$Length <- NULL
#Ctrl is your control or baseline sequences and the exp is the experimental sequences. 
names(CountTable)[1] <- "Ctrl1"
names(CountTable)[2] <- "Ctrl2"
names(CountTable)[3] <- "Ctrl3"
names(CountTable)[4] <- "EXP1"
names(CountTable)[5] <- "EXP2"
names(CountTable)[6] <- "EXP3"
samples <- data.frame(row.names=c("Ctrl1","Ctrl2","Ctrl3","EXP1","EXP2","EXP3"),condition=as.factor(c(rep("Ctrl",3), rep("EXP",3))))
Dataset <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~condition)


### Run DESEQ and print out results table, completed with gene annotations
DatasetProcessed <- DESeq(Dataset) # runs DESEQ
plotMA(DatasetProcessed, main="DESeq2", ylim=c(-5,5))

############## INDIVIDUAL CONTRAST ANALYSES ############## 

# Set Ctrl dataset as baseline
Dataset$condition <- relevel(Dataset$condition, "Ctrl")

res1 <- lfcShrink(DatasetProcessed, contrast=c("condition","CtrlEXP","Ctrl"))
baseMeanCtrl = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "Ctrl"])
baseMeanCtrlEXP = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "CtrlEXP"])
res1 = cbind(as.data.frame(res1), baseMeanCtrl, baseMeanCtrlEXP)
#we dont need to pull genes since featureCounts already labels them
#res1$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="SYMBOL", keytype="ENSEMBL", multiVals="first") # MAPS GENE IDs
#res1$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
write.csv(res1, "/Users/zachariasawaged/Documents/MS/bioinformatics/Assignments/Assignment_2_3/Assignment3/set_ctrl.csv", row.names=TRUE)

# GENERATE HISTOGRAM PLOTS OF pvalue distributions to confirm distribution is appropriate

hist(res1$pvalue, col = "lavender", main = "CtrlEXP vs Ctrl", xlab = "p-values")
hist(res1$pvalue[res1$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white") #this one states that genes with average count <1 should not be included in the graph

hist(res1$padj, col = "lavender", main = "CtrlEXP vs Ctrl", xlab = "padj-values")
hist(res1$padj[res1$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white") #this one states that genes with average count <1 should not be included in the graph


############## PLOTTING ############## 

### VOLCANO PLOTS
Dataset$condition <- relevel(Dataset$condition, "Ctrl")
res <- lfcShrink(DatasetProcessed, contrast=c("condition","EXP","Ctrl"))
with(res, plot(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, main="Volcano plot", xlim=c(-5,5)))
with(subset(res, padj>0.01), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="gray"))
with(subset(res, padj<0.01 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="red"))
with(subset(res, padj<0.01 & log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="blue"))
### ADD BELLS AND WHISTLES
pval = 0.01
abline(h = -log10(pval), col = "black", lty = 2, lwd=4)
mtext(paste("pval = 0.01", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 


### Plot heatmap - top20 variable genes
rld <- rlog(DatasetProcessed)
df <- as.data.frame(colData(DatasetProcessed))
topVarGenes <- head( order( rowMeans( assay(rld) ), decreasing=TRUE ), 20 )
mat <- assay(rld)[topVarGenes, ]
mat<-mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df) 


### Plot heatmap - selected genes (input file should contain list of Ensembl geneIDs, one per line)
mat <- read.table("LIST_OF_GENES.txt", header=TRUE, row.names=1)
row.names(mat) <- gene_list ##you miss this line
mat<-mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("external_gene_name","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
pheatmap(mat, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),border_color="NA")


### Plot heatmap - sorting by padj and log2fc
# PICK ALL GENES WITH pADJ < 0.05 AND THEN SUBSET FOR THOSE WITH Log2FC > 1 THEN PICK TOP 25 HITS
subset <- head((subset(res, res$log2FoldChange > 1 & res$padj < 0.00001)), n=25)  
sigGenes <- rownames(subset)
rows <- match(sigGenes, row.names(rld))
mat <- assay(rld)[rows,]
mat<-mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
df <- as.data.frame(colData(DatasetProcessed))
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df) 

# Exploratory plotting of count data using DESeq2 #

### LIBRARIES

library("ggplot2")

### Let's do some exploratory plotting. To show the effect of the transformation, we plot the first 
# sample against the second, first simply using the log2 function (after adding 1, to avoid taking the
# log of zero), and then using the rlog-transformed values. For the log2 method, we need estimate size 
# factors to account for sequencing depth (this is done automatically for the rlog method).
rlog <- rlog(Dataset)
par( mfrow = c( 1, 2 ) )
Dataset <- estimateSizeFactors(Dataset)
plot(log2( 1 + counts(Dataset, normalized=TRUE)[ , 1:2] ), pch=16, cex=0.3)
plot(assay(rlog)[ , 1:2], pch=16, cex=0.3)

###A useful first step in an RNA-Seq analysis is often to assess overall similarity between samples: 
# Which samples are similar to each other, which are different? Does this fit to the expectation from 
# the experiment's design?
sampleDists <- dist( t( assay(rlog) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

### Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). 
# In this ordination method, the data points (i.e., here, the samples) are projected onto the 
# 2D plane such that they spread out in the two directions which explain most of the differences 
# in the data. The x-axis is the direction (or principal component) which separates the data points
# the most. The amount of the total variance which is contained in the direction is printed in the 
# axis label.

# A useful exploration is to examine PCA plots generated from simple log transformted data (ntd1),
# rlog transformed data (rld1, our favorite), and VST transformed data (vst1, generally  used 
# for datasets of 20+ samples)

ntd1 <- normTransform(Dataset) #log2(n+1) transformation, the data that are transformed are read counts normalized for transcript average length and library size
rld1 <- rlog(Dataset, blind=FALSE) #the data that are transformed are read counts normalized for transcript average length and library size 
vsd1 <- varianceStabilizingTransformation(Dataset, blind=FALSE) #the data that are transformed are read counts normalized for transcript average length and library size 

# plot PCA for ntd1
pcaData <- plotPCA(ntd1, returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

# plot PCA for vsd1
pcaData <- plotPCA(vsd1, returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

# plot PCA for rld1
pcaDatar <- plotPCA(rld1, returnData=TRUE)
percentVarr <- round(100*attr(pcaDatar, "percentVar"))
ggplot(pcaDatar, aes(PC1, PC2, color=condition)) + geom_point(size =3) + xlab(paste0("PC1: ", percentVarr[1], "% variance")) + ylab(paste0("PC2: ", percentVarr[2], "% variance")) + coord_fixed()



