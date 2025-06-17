#########Downstream processing of count data in R

###reading data into R

CSrab18T1<-read.table("htseq_CSrab-18C-T1.txt", header = FALSE, col.names = c("Gene", "CSrab18T1"))
CSrab18T2<-read.table("htseq_CSrab-18C-T2.txt", header = FALSE, col.names = c("Gene", "CSrab18T2"))
CSrab18T3<-read.table("htseq_CSrab-18C-T3.txt", header = FALSE, col.names = c("Gene", "CSrab18T3"))
CSrab29T1<-read.table("htseq_CSrab-29C-T1.txt", header = FALSE, col.names = c("Gene", "CSrab29T1"))
CSrab29T2<-read.table("htseq_CSrab-29C-T2.txt", header = FALSE, col.names = c("Gene", "CSrab29T2"))
CSrab29T3<-read.table("htseq_CSrab-29C-T3.txt", header = FALSE, col.names = c("Gene", "CSrab29T3"))
IsoCS18T1<-read.table("htseq_IsoCS-18C-T1.txt", header = FALSE, col.names = c("Gene", "IsoCS18T1"))
IsoCS18T2<-read.table("htseq_IsoCS-18C-T2.txt", header = FALSE, col.names = c("Gene", "IsoCS18T2"))
IsoCS18T3<-read.table("htseq_IsoCS-18C-T3.txt", header = FALSE, col.names = c("Gene", "IsoCS18T3"))
IsoCS29T1<-read.table("htseq_IsoCS-29C-T1.txt", header = FALSE, col.names = c("Gene", "IsoCS29T1"))
IsoCS29T2<-read.table("htseq_IsoCS-29C-T2.txt", header = FALSE, col.names = c("Gene", "IsoCS29T2"))
IsoCS29T3<-read.table("htseq_IsoCS-29C-T3.txt", header = FALSE, col.names = c("Gene", "IsoCS29T3"))

###merging data into a single table

all.htseq__testes<- merge(CSrab18T1, CSrab18T2, by = "Gene",all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, CSrab18T3, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, CSrab29T1, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, CSrab29T2, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, CSrab29T3, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS18T1, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS18T2, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS18T3, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS29T1, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS29T2, by = "Gene", all = TRUE)
all.htseq__testes<- merge(all.htseq__testes, IsoCS29T3, by = "Gene", all = TRUE)

library(DESeq2)

###removing first five rows that contains data not being considered for analysis

all.htseq__testes<- all.htseq__testes[-1:-5, ]

### Changing rownames

rownames(all.htseq__testes)<- all.htseq__testes$Gene

### making a genotype array

Genotype<- c(rep("IIIm", 6), rep("Ym", 6))

### making a temperature array

Temperature<- c(rep("18C", 3), rep("29C", 3), rep("18C",3), rep("29C", 3))

### making a table containing samples, genotype and temperature

exp.design<- data.frame(Genotype, Temperature, rownames = names(all.htseq__testes[-c(1)]))

### DESeq analysis

dds<-DESeqDataSetFromMatrix(as.matrix(all.htseq__testes[,2:13]), colData = exp.design, design = ~ Genotype + Temperature + Genotype: Temperature)
dds<- DESeq(dds)
resultsNames(dds)

### Interaction  term (GxT): Temperature effect across genotypes

res<- results(dds, name = "GenotypeYm.Temperature29C")
res.df <- as.data.frame(res)
res_GT<- subset(res.df, padj<0.05)  ##Significantly differentially expressed genes
write.csv(res_GT, "interaction_testes.csv")

### Contrasting genotypes at 18C

geno1<- results(dds, contrast = c("Genotype", "Ym", "IIIm"))
geno1.df<- data.frame(geno1)
geno_18C<- subset(geno1.df, padj<0.05)  ##Significantly differentiated expressed genes 

### Contrasting genotypes at 29C

geno2<- results(dds, list(c("Genotype_Ym_vs_IIIm", "GenotypeYm.Temperature29C")))
geno2<- data.frame(geno2)
geno_29C<- data.frame(subset(geno2, padj<0.05))   ##Significantly differentially expressed genes

### Effect of temperature on genotype IIIM

res1<- results(dds, contrast = c("Temperature", "29C", "18C"))
res1.df <- data.frame(res1)
temp_IIIM<- subset(res1.df, padj<0.05)   ##Significantly differentiated expressed genes

### Effect of temperature on genotype YM

res2<- results(dds, contrast =  list(c("Temperature_29C_vs_18C", "GenotypeYm.Temperature29C")))
res2.df<- as.data.frame(res2)
temp_YM<- subset(res2.df, padj<0.05)    ##Significantly differentiated expressed genes


### Plotting PCA

rld<- rlog(dds)
plotPCA(rld, intgroup= c("Temperature", "Genotype"), ntop=16540) ###Using all genes
plotPCA(rld, intgroup= c("Temperature", "Genotype"), ntop=500)  ### Using top 500 genes

### Heatmap of sample to sample distance
##calculate the distance

rld<- rlog(dds)
Dis<- dist(t(assay(rld)))

##heatmap
library("RColorBrewer")

DistMatrix <- as.matrix(Dis)
rownames(DistMatrix) <- paste(rld$Temperature, rld$Genotype, sep="-")
colnames(DistMatrix) <-  paste(rld$Temperature, rld$Genotype, sep="-")
colors <- colorRampPalette((brewer.pal(9, "Blues")) )(255)      # 0-white and 100-blue , use rev function to reverse it
library("pheatmap")
pheatmap(DistMatrix,
         clustering_distance_rows=Dis,
         clustering_distance_cols=Dis,
         col=colors, treeheight_col = 0)      # treeheight_col=0 removes the dendograms at the top


### For MA plot with lfc shrinkage
## as type="apeglm" can not be used in contrasts, I will be using type ="ashr"
### GxT effect
resLFC <- lfcShrink(dds, coef="GenotypeYm.Temperature29C", type="ashr")
plotMA(resLFC, ylim=c(-5,5), alpha=0.05)

### Genotype effect at 18C
geno1_LFC<- lfcShrink(dds, coef = "Genotype_Ym_vs_IIIm", type = "ashr")
plotMA(geno1_LFC, ylim=c(-5,5), alpha=0.05)

### Genotype effect at 29C
geno2_LFC<- lfcShrink(dds, contrast= list(c("Genotype_Ym_vs_IIIm", "GenotypeYm.Temperature29C")),type = "ashr")
plotMA(geno2_LFC, ylim=c(-5,5), alpha=0.05)

### Temperature effect in IIIM
res2_LFC<- lfcShrink(dds, coef ="Temperature_29C_vs_18C", type = "ashr")
plotMA(res2_LFC, ylim=c(-5,5), alpha=0.05)

### Temperature effect in YM
res3_LFC<- lfcShrink(dds, contrast =  list(c("Temperature_29C_vs_18C", "GenotypeYm.Temperature29C")), type="ashr")
plotMA(res3_LFC, ylim=c(-5,5), alpha=0.05)



