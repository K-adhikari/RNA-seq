#### #####Downstream processing of count data in R

#reading data into R
CSrab18H1<-read.table("htseq_CSrab-18C-H1.txt", header = FALSE, col.names = c("Gene", "CSrab18H1"))
CSrab18H2<-read.table("htseq_CSrab-18C-H2.txt", header = FALSE, col.names = c("Gene", "CSrab18H2"))
CSrab18H3<-read.table("htseq_CSrab-18C-H3.txt", header = FALSE, col.names = c("Gene", "CSrab18H3"))
CSrab29H1<-read.table("htseq_CSrab-29C-H1.txt", header = FALSE, col.names = c("Gene", "CSrab29H1"))
CSrab29H2<-read.table("htseq_CSrab-29C-H2.txt", header = FALSE, col.names = c("Gene", "CSrab29H2"))
CSrab29H3<-read.table("htseq_CSrab-29C-H3.txt", header = FALSE, col.names = c("Gene", "CSrab29H3"))
IsoCS18H1<-read.table("htseq_IsoCS-18C-H1.txt", header = FALSE, col.names = c("Gene", "IsoCS18H1"))
IsoCS18H2<-read.table("htseq_IsoCS-18C-H2.txt", header = FALSE, col.names = c("Gene", "IsoCS18H2"))
IsoCS18H3<-read.table("htseq_IsoCS-18C-H3.txt", header = FALSE, col.names = c("Gene", "IsoCS18H3"))
IsoCS29H1<-read.table("htseq_IsoCS-29C-H1.txt", header = FALSE, col.names = c("Gene", "IsoCS29H1"))
IsoCS29H2<-read.table("htseq_IsoCS-29C-H2.txt", header = FALSE, col.names = c("Gene", "IsoCS29H2"))
IsoCS29H3<-read.table("htseq_IsoCS-29C-H3.txt", header = FALSE, col.names = c("Gene", "IsoCS29H3"))

#merging data into a single table
all.htseq__head<- merge(CSrab18H1, CSrab18H2, by = "Gene",all = TRUE)
all.htseq__head<- merge(all.htseq__head, CSrab18H3, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, CSrab29H1, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, CSrab29H2, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, CSrab29H3, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS18H1, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS18H2, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS18H3, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS29H1, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS29H2, by = "Gene", all = TRUE)
all.htseq__head<- merge(all.htseq__head, IsoCS29H3, by = "Gene", all = TRUE)

library(DESeq2)

## removing first five rows that contains data not being considered for analysis

all.htseq__head<- all.htseq__head[-1:-5, ]

## Changing rownames
rownames(all.htseq__head)<- all.htseq__head$Gene

## making a genotype array
Genotype<- c(rep("IIIm", 6), rep("Ym", 6))

## making a temperature array
Temperature<- c(rep("18C", 3), rep("29C", 3), rep("18C",3), rep("29C", 3))

## making a table containing samples, genotype and temperature
exp.design<- data.frame(Genotype, Temperature, rownames = names(all.htseq__head[-c(1)]))

## DESeq analysis
dds<-DESeqDataSetFromMatrix(as.matrix(all.htseq__head[,2:13]), colData = exp.design, design = ~ Genotype + Temperature + Genotype: Temperature)
dds<- DESeq(dds)
resultsNames(dds)

## Interaction  term (GxT): Temperature effect across genotypes 
res1<- data.frame(results(dds, name="GenotypeYm.Temperature29C"))
res_GT<- subset(res1, padj<0.05)
write.csv(res_GT, "interaction_head.csv")

## Contrasting genotypes at 18C
geno1<- results(dds, contrast = c("Genotype", "Ym", "IIIm"))
geno1.df <- data.frame(geno1)
geno_18C<- subset(geno1.df, padj<0.05)

## Contrasting genotypes at 29C
geno2<- results(dds, list(c("Genotype_Ym_vs_IIIm", "GenotypeYm.Temperature29C")))
geno2.df<- data.frame(geno2)
geno_29C<- subset(geno2.df, padj<0.05)

## Effect of temperature on genotype IIIm
res2<- results(dds, contrast = c("Temperature", "29C", "18C"))
res2.df <- as.data.frame(res2)
temp_IIIM<- subset(res2.df, padj<0.05)

#Effect of temperature on genotype Ym
res3<- results(dds, contrast =  list(c("Temperature_29C_vs_18C", "GenotypeYm.Temperature29C")))
res3.df <- as.data.frame(res3)
temp_YM<- subset(res3.df, padj<0.05)

## Plotting PCA
rld<- rlog(dds)
plotPCA(rld, intgroup= c("Temperature", "Genotype"), ntop=16540) ###Using all genes
plotPCA(rld, intgroup= c("Temperature", "Genotype"), ntop=500)  ### Using top 500 genes


## Heatmap of sample to sample distance

#calculate the distance
rld<- rlog(dds) 
Dis<- dist(t(assay(rld)))

#heatmap

library("RColorBrewer")

sampleDistMatrix <- as.matrix(Dis)
rownames(sampleDistMatrix) <- paste(rld$Temperature, rld$Genotype, sep="-")
colnames(sampleDistMatrix) <- paste(rld$Temperature, rld$Genotype, sep="-")
colors <- colorRampPalette((brewer.pal(9, "Blues")))(255)  # 0-white and 100-blue , use rev function to reverse it
library("pheatmap")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=Dis,
         clustering_distance_cols=Dis,
         col=colors, treeheight_col = 0) # treeheight_col=0 removes the dendograms at the top.


### For MA plot with lfc shrinkage
## as type="apeglm" can not be used in contrasts, I will be using type ="ashr"
## For GxT effect

resLFC <- lfcShrink(dds, coef="GenotypeYm.Temperature29C", type="ashr")
plotMA(resLFC, ylim=c(-5,5), alpha=0.05)

### For Genotype effect at 18C
geno1_LFC<- lfcShrink(dds, coef = "Genotype_Ym_vs_IIIm", type = "ashr")
plotMA(geno1_LFC, ylim=c(-5,5), alpha=0.05)

#### For Genotype effect at 29C
geno2_LFC<- lfcShrink(dds, contrast= list(c("Genotype_Ym_vs_IIIm", "GenotypeYm.Temperature29C")),type = "ashr")
plotMA(geno2_LFC, ylim=c(-5,5), alpha=0.05)

#### For Temperature effect in IIIM
res2_LFC<- lfcShrink(dds, coef ="Temperature_29C_vs_18C", type = "ashr")
plotMA(res2_LFC, ylim=c(-5,5), alpha=0.05)

##### For Temperature effect in YM
res3_LFC<- lfcShrink(dds, contrast =  list(c("Temperature_29C_vs_18C", "GenotypeYm.Temperature29C")), type="ashr")
plotMA(res3_LFC, ylim=c(-5,5), alpha=0.05)
