### Downstream processing of filtered variants
library(vcfR)
vcf_testis<- read.vcfR("Testis_JointGenotype_Filtered_Pass.vcf", verbose= FALSE)
gt<- extract.gt(vcf_testis)

#checking heterozygous positions
hets<- is_het(gt)
#changing into a dataframe
hets<- as.data.frame(hets)
#Giving name to first column
hets<- setNames(cbind(rownames(hets), hets, row.names=NULL), c("Replicon Accession", "CSrab18Ctestis", "CSrab29Ctestis", "IsoCS18Ctestis", "IsoCS29Ctestis"))

# heterozygous common for IIIm at 18 and 29C and not in Ym
hets_IIIm_testis<- subset(hets, hets$CSrab18Ctestis==TRUE & hets$CSrab29Ctestis==TRUE & hets$IsoCS18Ctestis==FALSE & hets$IsoCS29Ctestis==FALSE)
hets_IIIm_testis<- hets_IIIm_testis[-c(2:5)]

#Extracting alleles
A<- extract.gt(vcf_testis, return.alleles=TRUE)
#Filtering out indels
B<- extract.indels(vcf_testis)
#alleles without indels
B1<- extract.gt(B, return.alleles = TRUE)

#Creating a data frame from all SNPs
B2<- as.data.frame(B1, header=TRUE)

#Giving name to first column
B2<- setNames(cbind(rownames(B2), B2, row.names=NULL), c("Replicon Accession", "CSrab18Ctestis", "CSrab29Ctestis", "IsoCS18Ctestis", "IsoCS29Ctestis"))

#merge data frame B2 and hets_IIIm
IIIm_testis<- merge(hets_IIIm_testis, B2, by="Replicon Accession")

#separating position from scaffold accession
library(tidyr)
IIIm_testis<- separate(IIIm_testis, `Replicon Accession`, into=c("NW", "Replicon.Accession", "Position"), sep="_")
IIIm_testis$'Replicon Accession'<- paste(IIIm_testis$NW,IIIm_testis$`Replicon.Accession`, sep="_")
IIIm_testis<- IIIm_testis[-c(1:2)]
IIIm_testis<- IIIm_testis[c(6,1,2,3,4,5)]
IIIm_testis$Position<- as.numeric(IIIm_testis$Position)
IIIm_testis$`Replicon Accession`<- as.factor(IIIm_testis$`Replicon Accession`)
library(dplyr)
IIIm_testis<- IIIm_testis%>%arrange(`Replicon Accession`, Position)
IIIm_testis$order<- c(1:26529)


#Chromosome 1
chr1_scaff<- read.csv("chromosome1_scaffold.csv", header = TRUE)
chr1_scaff<- chr1_scaff[c(-1)]
colnames(chr1_scaff)[1]<- "Replicon Accession"
chr1_IIIm_testis<- merge(chr1_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chr1_IIIm_testis<- chr1_IIIm_testis%>%arrange(order)


#Chromosome 2
chr2_scaff<- read.csv("chromosome2_scaffold.csv", header = TRUE)
chr2_scaff<- chr2_scaff[c(-1)]
colnames(chr2_scaff)[1]<- "Replicon Accession"
chr2_IIIm_testis<- merge(chr2_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chr2_IIIm_testis<- chr2_IIIm_testis%>%arrange(order)


#Chromosome 3
chr3_scaff<- read.csv("chromosome3_scaffold.csv", header = TRUE)
chr3_scaff<- chr3_scaff[c(-1)]
colnames(chr3_scaff)[1]<- "Replicon Accession"
chr3_IIIm_testis<- merge(chr3_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chr3_IIIm_testis<- chr3_IIIm_testis%>%arrange(order)


#Chromosome 4
chr4_scaff<- read.csv("chromosome4_scaffold.csv", header = TRUE)
chr4_scaff<- chr4_scaff[c(-1)]
colnames(chr4_scaff)[1]<- "Replicon Accession"
chr4_IIIm_testis<- merge(chr4_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chr4_IIIm_testis<- chr4_IIIm_testis%>%arrange(order)


#Chromosome 5
chr5_scaff<- read.csv("chromosome5_scaffold.csv", header = TRUE)
chr5_scaff<- chr5_scaff[c(-1)]
colnames(chr5_scaff)[1]<- "Replicon Accession"
chr5_IIIm_testis<- merge(chr5_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chr5_IIIm_testis<- chr5_IIIm_testis%>%arrange(order)


#Chromosome X
chrX_scaff<- read.csv("chromosomeX_scaffold.csv", header = TRUE)
chrX_scaff<- chrX_scaff[c(-1)]
colnames(chrX_scaff)[1]<- "Replicon Accession"
chrX_IIIm_testis<- merge(chrX_scaff, IIIm_testis, by="Replicon Accession")
library(dplyr)
chrX_IIIm_testis<- chrX_IIIm_testis%>%arrange(order)


write.csv(chr1_IIIm_testis, "1stchromosome_hetsIIIm_testis.csv")
write.csv(chr2_IIIm_testis, "2ndchromosome_hetsIIIm_testis.csv")
write.csv(chr3_IIIm_testis, "3rdchromosome_hetsIIIm_testis.csv")
write.csv(chr4_IIIm_testis, "4thchromosome_hetsIIIm_testis.csv")
write.csv(chr5_IIIm_testis, "5thchromosome_hetsIIIm_testis.csv")
write.csv(chrX_IIIm_testis, "Xchromosome_hetsIIIm_testis.csv")



####### Geetting gene info for the scaffolds. 
# This portion of code should be run on a high performance computing cluster as this takes a long time to do 

## 1st Chromosome
C1<- read.csv("1stchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
C1$'Replicon.Accession'<- as.character(C1$'Replicon.Accession')
C1$Replicon.Accession<- gsub("\\..*","",C1$Replicon.Accession)
C1$Position<- as.numeric(C1$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C1 for gene names
C1$gene = NA
#Getting information from source to file with allele information
fun<- function(x){
  for (i in 1:2060){
    for(k in 1:17508){
      if(C1[i,2]==source[k,7]){
        if(C1[i,3]>=source[k,8] & C1[i,3]<=source[k,9]){
          C1[i,8]=source[k,5]
          C1
        }
      }
    }
  }
  return (C1)
}
C1_1<- fun()
write.csv(C1,"IIImhets_1stchromosome_with_gene_info.csv")


### 2nd chromosome

C2<- read.csv("2ndchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
C2$'Replicon.Accession'<- as.character(C2$'Replicon.Accession')
C2$Replicon.Accession<- gsub("\\..*","",C2$Replicon.Accession)
C2$Position<- as.numeric(B2$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C2 for gene names
C2$gene = NA
#Getting information from source to file with allele information
fun1<- function(x){
  for (i in 1:3254){
    for(k in 1:17508){
      if(C2[i,2]==source[k,7]){
        if(C2[i,3]>=source[k,8] & C2[i,3]<=source[k,9]){
          C2[i,9]=source[k,5]
          C2
        }
      }
    }
  }
  return (C2)
}
C2_1<- fun1()
write.csv(C2_1,"IIImhets_2ndchromosome_with_gene_info.csv")


## 3rd chromosome

C3<- read.csv("3rdchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
C3$'Replicon.Accession'<- as.character(C3$'Replicon.Accession')
C3$Replicon.Accession<- gsub("\\..*","",C3$Replicon.Accession)
C3$Position<- as.numeric(C3$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C3 for gene names
C3$gene = NA
#Getting information from source to file with allele information
fun2<- function(x){
  for (i in 1:13773){
    for(k in 1:17508){
      if(C3[i,2]==source[k,7]){
        if(C3[i,3]>=source[k,8] & C3[i,3]<=source[k,9]){
          C3[i,9]=source[k,5]
          C3
        }
      }
    }
  }
  return (C3)
}
C3_1<- fun2()
write.csv(C3_1,"IIImhets_3rdchromosome_with_gene_info.csv")

## 4th chromosome

C4<- read.csv("4thchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
C4$'Replicon.Accession'<- as.character(C4$'Replicon.Accession')
C4$Replicon.Accession<- gsub("\\..*","",C4$Replicon.Accession)
C4$Position<- as.numeric(C4$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C4 for gene names
C4$gene = NA
#Getting information from source to file with allele information
fun3<- function(x){
  for (i in 1:1863){
    for(k in 1:17508){
      if(C4[i,2]==source[k,7]){
        if(C4[i,3]>=source[k,8] & C4[i,3]<=source[k,9]){
          C4[i,9]=source[k,5]
          C4
        }
      }
    }
  }
  return (C4)
}
C4_1<- fun3()
write.csv(C4_1,"IIImhets_4thchromosome_with_gene_info.csv")

## 5th chromosome

C5<- read.csv("5thchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
C5$'Replicon.Accession'<- as.character(C5$'Replicon.Accession')
C5$Replicon.Accession<- gsub("\\..*","",C5$Replicon.Accession)
C5$Position<- as.numeric(C5$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C5 for gene names
C5$gene = NA
#Getting information from source to file with allele information
fun4<- function(x){
  for (i in 1:2586){
    for(k in 1:17508){
      if(C5[i,2]==source[k,7]){
        if(C5[i,3]>=source[k,8] & C5[i,3]<=source[k,9]){
          C5[i,9]=source[k,5]
          C5
        }
      }
    }
  }
  return (C5)
}
C5_1<- fun4()
write.csv(C5_1,"IIImhets_5thchromosome_with_gene_info.csv")


## X chromosome

CX<- read.csv("Xchromosome_hetsIIIm_testis.csv", header=TRUE, sep=",")
CX$'Replicon.Accession'<- as.character(CX$'Replicon.Accession')
CX$Replicon.Accession<- gsub("\\..*","",CX$Replicon.Accession)
CX$Position<- as.numeric(CX$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for CX for gene names
CX$gene = NA
#Getting information from source to file with allele information
fun4<- function(x){
  for (i in 1:77){
    for(k in 1:17508){
      if(CX[i,2]==source[k,7]){
        if(CX[i,3]>=source[k,8] & CX[i,3]<=source[k,9]){
          CX[i,9]=source[k,5]
          CX
        }
      }
    }
  }
  return (CX)
}
CX_1<- fun5()
write.csv(CX_1,"IIImhets_Xchromosome_with_gene_info.csv")




######## Extracting allele depth from vcf files
ad<- extract.gt(vcf, element = "AD")
ad<- as.data.frame(ad)
ad<- setNames(cbind(rownames(ad), ad, row.names=NULL), c("scaff", "CSrab18Ctestis", "CSrab29Ctestis", "IsoCS18Ctestis", "IsoCS29Ctestis"))
ad$scaff<- as.character(ad$scaff)
library(dplyr)
library(tidyr)
ad<- separate(ad, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep=",")
ad<- separate(ad, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep=",")
ad<- separate(ad, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep=",")
ad<- separate(ad, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep=",")
ad$CSrab18Ctestis.1<- as.numeric(ad$CSrab18Ctestis.1)
ad$CSrab18Ctestis.2<- as.numeric(ad$CSrab18Ctestis.2)
ad$CSrab29Ctestis.1<- as.numeric(ad$CSrab29Ctestis.1)
ad$CSrab29Ctestis.2<- as.numeric(ad$CSrab29Ctestis.2)
ad$IsoCS18Ctestis.1<- as.numeric(ad$IsoCS18Ctestis.1)
ad$IsoCS18Ctestis.2<- as.numeric(ad$IsoCS18Ctestis.2)
ad$IsoCS29Ctestis.1<- as.numeric(ad$IsoCS29Ctestis.1)
ad$IsoCS29Ctestis.2<- as.numeric(ad$IsoCS29Ctestis.2)

#normalizing allele depths
ad$CSrab18Ctestis.1<- ad$CSrab18Ctestis.1*1000000/61366316
ad$CSrab18Ctestis.2<- ad$CSrab18Ctestis.2*1000000/61366316
ad$CSrab29Ctestis.1<- ad$CSrab29Ctestis.1*1000000/98504222
ad$CSrab29Ctestis.2<- ad$CSrab29Ctestis.2*1000000/98504222
ad$IsoCS18Ctestis.1<- ad$IsoCS18Ctestis.1*1000000/95559868
ad$IsoCS18Ctestis.2<- ad$IsoCS18Ctestis.2*1000000/95559868
ad$IsoCS29Ctestis.1<- ad$IsoCS29Ctestis.1*1000000/117629242
ad$IsoCS29Ctestis.2<- ad$IsoCS29Ctestis.2*1000000/117629242

#Putting YM together
ad1<- within(ad, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ad1<- within(ad1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ad1<-ad1[, c(1:5,10,11)]

library(reshape2)
ad<- melt(ad, id="scaff")

#For new[1]dataframe
ad1<- melt(ad1, id="scaff")


####For genes on the first chromosome


read_IIIm1<- read.csv("IIImhets_1stchromosome_with_gene_info.csv", header=TRUE, sep = ",")

###Reading genes with significant GxT interaction from DESeq2 analysis
interaction<- read.csv("interaction_testes.csv", header=TRUE, sep=",")

####In order to limit the analysis to DE genes
ASE_IIImhets<- merge(interaction, read_IIIm1, by="gene")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:280,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:210,]
#ASE_IIImhets1<- ASE_IIImhets1[1:444,]

New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                      "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis",
                                      "IsoCS18Ctestis.1"="IsoCS18Ctestis", "IsoCS18Ctestis.2"="IsoCS18Ctestis",
                                      "IsoCS29Ctestis.1"="IsoCS29Ctestis", "IsoCS29Ctestis.2"="IsoCS29Ctestis")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,280,by=8))
    for(k in seq(1,210, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Ctestis" & New1[k+5,2]=="IsoCS29Ctestis"){
            New1[k+4,5]= New[i+4,5] + New[i+5,5]
            New1[k+5,5]= New[i+6,5] + New[i+7,5]
            New1
          }
        }
      }
    }
  return(New1)
}

New2<- data.frame(fun())
New2$Normalized_count<- as.numeric(New2$Normalized_count)
New2<- transform(New2, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                        "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis")))

write.csv(New2, "Testis_1st_chromosome_heterozygous_data.csv")


#### For genes on the second chromosome

read_IIIm2<- read.csv("IIImhets_2ndchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIIm2, by="gene")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:176,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:132,]


New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                      "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis",
                                      "IsoCS18Ctestis.1"="IsoCS18Ctestis", "IsoCS18Ctestis.2"="IsoCS18Ctestis",
                                      "IsoCS29Ctestis.1"="IsoCS29Ctestis", "IsoCS29Ctestis.2"="IsoCS29Ctestis")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,176,by=8))
    for(k in seq(1,132, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Ctestis" & New1[k+5,2]=="IsoCS29Ctestis"){
            New1[k+4,5]= New[i+4,5] + New[i+5,5]
            New1[k+5,5]= New[i+6,5] + New[i+7,5]
            New1
          }
        }
      }
    }
  return(New1)
}

New2<- data.frame(fun())
New2$Normalized_count<- as.numeric(New2$Normalized_count)
New2<- transform(New2, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                        "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis")))

write.csv(New2, "Testis_2nd_chromosome_heterozygous_data.csv")



#### For genes on the third chromosome

read_IIIm3<- read.csv("IIImhets_3rdchromosome_with_gene_info.csv", header=TRUE, sep = ",")


ASE_IIImhets<- merge(interaction, read_IIIm3, by="gene")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:760,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:570,]


New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                      "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis",
                                      "IsoCS18Ctestis.1"="IsoCS18Ctestis", "IsoCS18Ctestis.2"="IsoCS18Ctestis",
                                      "IsoCS29Ctestis.1"="IsoCS29Ctestis", "IsoCS29Ctestis.2"="IsoCS29Ctestis")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,760,by=8))
    for(k in seq(1,570, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Ctestis" & New1[k+5,2]=="IsoCS29Ctestis"){
            New1[k+4,5]= New[i+4,5] + New[i+5,5]
            New1[k+5,5]= New[i+6,5] + New[i+7,5]
            New1
          }
        }
      }
    }
  return(New1)
}

New2<- data.frame(fun())
New2$Normalized_count<- as.numeric(New2$Normalized_count)
New2<- transform(New2, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                        "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis")))

write.csv(New2, "Testis_3rd_chromosome_heterozygous_data.csv")


#### For fourth chromosome genes

read_IIIm4<- read.csv("IIImhets_4thchromosome_with_gene_info.csv", header=TRUE, sep = ",")


ASE_IIImhets<- merge(interaction, read_IIIm4, by="gene")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:136,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:102,]


New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                      "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis",
                                      "IsoCS18Ctestis.1"="IsoCS18Ctestis", "IsoCS18Ctestis.2"="IsoCS18Ctestis",
                                      "IsoCS29Ctestis.1"="IsoCS29Ctestis", "IsoCS29Ctestis.2"="IsoCS29Ctestis")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,136,by=8))
    for(k in seq(1,102, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Ctestis" & New1[k+5,2]=="IsoCS29Ctestis"){
            New1[k+4,5]= New[i+4,5] + New[i+5,5]
            New1[k+5,5]= New[i+6,5] + New[i+7,5]
            New1
          }
        }
      }
    }
  return(New1)
}

New2<- data.frame(fun())
New2$Normalized_count<- as.numeric(New2$Normalized_count)
New2<- transform(New2, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                        "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis")))

write.csv(New2, "Testis_4th_chromosome_heterozygous_data.csv")

  
#### For fifth chromosome genes


read_IIIm5<- read.csv("IIImhets_5thchromosome_with_gene_info.csv", header=TRUE, sep = ",")


ASE_IIImhets<- merge(interaction, read_IIIm5, by="gene")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Ctestis, into = c("CSrab18Ctestis.1", "CSrab18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Ctestis, into = c("CSrab29Ctestis.1", "CSrab29Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Ctestis, into = c("IsoCS18Ctestis.1", "IsoCS18Ctestis.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Ctestis, into = c("IsoCS29Ctestis.1", "IsoCS29Ctestis.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Ctestis<- paste(IsoCS18Ctestis.1, IsoCS18Ctestis.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Ctestis<- paste(IsoCS29Ctestis.1, IsoCS29Ctestis.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:296,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:222,]


New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                      "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis",
                                      "IsoCS18Ctestis.1"="IsoCS18Ctestis", "IsoCS18Ctestis.2"="IsoCS18Ctestis",
                                      "IsoCS29Ctestis.1"="IsoCS29Ctestis", "IsoCS29Ctestis.2"="IsoCS29Ctestis")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,296,by=8))
    for(k in seq(1,222, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Ctestis" & New1[k+5,2]=="IsoCS29Ctestis"){
            New1[k+4,5]= New[i+4,5] + New[i+5,5]
            New1[k+5,5]= New[i+6,5] + New[i+7,5]
            New1
          }
        }
      }
    }
  return(New1)
}

New2<- data.frame(fun())
New2$Normalized_count<- as.numeric(New2$Normalized_count)
New2<- transform(New2, GT=revalue(GT, c("CSrab18Ctestis.1"="CSrab18Ctestis","CSrab18Ctestis.2"="CSrab18Ctestis",
                                        "CSrab29Ctestis.1"="CSrab29Ctesis", "CSrab29Ctestis.2"="CSrab29Ctestis")))

write.csv(New2, "Testis_5th_chromosome_heterozygous_data.csv")



#### For X chromosome genes


read_IIImX<- read.csv("IIImhets_Xchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIImX, by="gene") #No common genes were found so no downstream analysis required



####################### IIIIM-III calculations
#### Chromosome 1


Total<- read.csv("Testis_1st_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total<- Total[-c(163:168),]
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),34)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),34)
Total$order<- 1:204
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:135){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:68){
        if (df1[k,2]==a & df1 [k,3]==b) {
          df1[k,4]= T3[i+1,6]- T3[i,6]
          df1
          
        }
      }
    }
  }
  return (df1)
}

df2<- data.frame(funct())

#To calculate average within a gene
df3<- df2[!duplicated(df2[c(1,3)]),]
df3<- df3[,c(1,3)]

df4<- data.frame(unique(df3$gene))
colnames(df4)[1]<- "gene"


new<- function(x){
  for (i in 1:4) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:8) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:8){
      if (df3[k,1]==a & df3[k,2]=="29C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }
    }}
  return(df3)
}

df6<- data.frame(new())
colnames(df6)[3]<- "difference"
colnames(df6)[4]<- "SE"

write.csv(df6, "Testis_1st_chromosome_gene_details")


#####Chromosome 2

Total<- read.csv("Testis_2nd_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total<- Total[-c(25:36,85:90),]
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),19)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),19)
Total$order<- 1:114
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:75){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:38){
        if (df1[k,2]==a & df1 [k,3]==b) {
          df1[k,4]= T3[i+1,6]- T3[i,6]
          df1
          
        }
      }
    }
  }
  return (df1)
}

df2<- data.frame(funct())

#To calculate average within a gene
df3<- df2[!duplicated(df2[c(1,3)]),]
df3<- df3[,c(1,3)]

df4<- data.frame(unique(df3$gene))
colnames(df4)[1]<- "gene"


new<- function(x){
  for (i in 1:10) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:20) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:20){
      if (df3[k,1]==a & df3[k,2]=="29C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }
    }}
  return(df3)
}

df6<- data.frame(new())
colnames(df6)[3]<- "difference"
colnames(df6)[4]<- "SE"

write.csv(df6, "Testis_2nd_chromosome_gene_details")


#####Chromosome 3

Total<- read.csv("Testis_3rd_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),95)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),95)
Total$order<- 1:570
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:379){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:190){
        if (df1[k,2]==a & df1 [k,3]==b) {
          df1[k,4]= T3[i+1,6]- T3[i,6]
          df1
          
        }
      }
    }
  }
  return (df1)
}

df2<- data.frame(funct())

#To calculate average within a gene
df3<- df2[!duplicated(df2[c(1,3)]),]
df3<- df3[,c(1,3)]

df4<- data.frame(unique(df3$gene))
colnames(df4)[1]<- "gene"


new<- function(x){
  for (i in 1:12) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:24) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:24){
      if (df3[k,1]==a & df3[k,2]=="29C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }
    }}
  return(df3)
}

df6<- data.frame(new())
colnames(df6)[3]<- "difference"
colnames(df6)[4]<- "SE"

write.csv(df6, "Testis_3rd_chromosome_gene_details")



####Chromosome 4

Total<- read.csv("Testis_4th_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),17)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),17)
Total$order<- 1:102
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:67){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:34){
        if (df1[k,2]==a & df1 [k,3]==b) {
          df1[k,4]= T3[i+1,6]- T3[i,6]
          df1
          
        }
      }
    }
  }
  return (df1)
}

df2<- data.frame(funct())

#To calculate average within a gene
df3<- df2[!duplicated(df2[c(1,3)]),]
df3<- df3[,c(1,3)]

df4<- data.frame(unique(df3$gene))
colnames(df4)[1]<- "gene"


new<- function(x){
  for (i in 1:4) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:8) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:8){
      if (df3[k,1]==a & df3[k,2]=="29C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }
    }}
  return(df3)
}

df6<- data.frame(new())
colnames(df6)[3]<- "difference"
colnames(df6)[4]<- "SE"

write.csv(df6, "Testis_4th_chromosome_gene_details")


##### Chromosome 5

Total<- read.csv("Testis_5th_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),37)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),37)
Total$order<- 1:222
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:147){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:74){
        if (df1[k,2]==a & df1 [k,3]==b) {
          df1[k,4]= T3[i+1,6]- T3[i,6]
          df1
          
        }
      }
    }
  }
  return (df1)
}

df2<- data.frame(funct())

#To calculate average within a gene
df3<- df2[!duplicated(df2[c(1,3)]),]
df3<- df3[,c(1,3)]

df4<- data.frame(unique(df3$gene))
colnames(df4)[1]<- "gene"


new<- function(x){
  for (i in 1:7) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:14) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:14){
      if (df3[k,1]==a & df3[k,2]=="29C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }
    }}
  return(df3)
}

df6<- data.frame(new())
colnames(df6)[3]<- "difference"
colnames(df6)[4]<- "SE"

write.csv(df6, "Testis_5th_chromosome_gene_details")
