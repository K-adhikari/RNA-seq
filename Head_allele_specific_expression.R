library(vcfR)
vcf_head<- read.vcfR("Head_JointGenotype_Filtered_Pass.vcf", verbose= FALSE)
gt<- extract.gt(vcf_head)

#checking heterozygous positions
hets<- is_het(gt)
#changing into a dataframe
hets<- as.data.frame(hets)
#Giving name to first column
hets<- setNames(cbind(rownames(hets), hets, row.names=NULL), c("Replicon Accession", "CSrab18Chead", "CSrab29Chead", "IsoCS18Chead", "IsoCS29Chead"))

# heterozygous common for IIIm at 18 and 29C and not in Ym
hets_IIIm_head<- subset(hets, hets$CSrab18Chead==TRUE & hets$CSrab29Chead==TRUE & hets$IsoCS18Chead==FALSE & hets$IsoCS29Chead==FALSE)
hets_IIIm_head<- hets_IIIm_head[-c(2:5)]

#Extracting alleles
A<- extract.gt(vcf_head, return.alleles=TRUE)
#Filtering out indels
B<- extract.indels(vcf_head)
#alleles without indels
B1<- extract.gt(B, return.alleles = TRUE)

#Creating a data frame from all SNPs
B2<- as.data.frame(B1, header=TRUE)

#Giving name to first column
B2<- setNames(cbind(rownames(B2), B2, row.names=NULL), c("Replicon Accession", "CSrab18Chead", "CSrab29Chead", "IsoCS18Chead", "IsoCS29Chead"))

#merge data frame B2 and hets_IIIm
IIIm_head<- merge(hets_IIIm_head, B2, by="Replicon Accession")

#separating position from scaffold accession
library(tidyr)
IIIm_head<- separate(IIIm_head, `Replicon Accession`, into=c("NW", "Replicon.Accession", "Position"), sep="_")
IIIm_head$'Replicon Accession'<- paste(IIIm_head$NW,IIIm_head$`Replicon.Accession`, sep="_")
IIIm_head<- IIIm_head[-c(1:2)]
IIIm_head<- IIIm_head[c(6,1,2,3,4,5)]
IIIm_head$Position<- as.numeric(IIIm_head$Position)
IIIm_head$`Replicon Accession`<- as.factor(IIIm_head$`Replicon Accession`)
library(dplyr)
IIIm_head<- IIIm_head%>%arrange(`Replicon Accession`, Position)
IIIm_head$order<- c(1:32061)


#Chromosome 1
chr1_scaff<- read.csv("chromosome1_scaffold.csv", header = TRUE)
chr1_scaff<- chr1_scaff[c(-1)]
colnames(chr1_scaff)[1]<- "Replicon Accession"
chr1_IIIm_head<- merge(chr1_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chr1_IIIm_head<- chr1_IIIm_head%>%arrange(order)


#Chromosome 2
chr2_scaff<- read.csv("chromosome2_scaffold.csv", header = TRUE)
chr2_scaff<- chr2_scaff[c(-1)]
colnames(chr2_scaff)[1]<- "Replicon Accession"
chr2_IIIm_head<- merge(chr2_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chr2_IIIm_head<- chr2_IIIm_head%>%arrange(order)


#Chromosome 3
chr3_scaff<- read.csv("chromosome3_scaffold.csv", header = TRUE)
chr3_scaff<- chr3_scaff[c(-1)]
colnames(chr3_scaff)[1]<- "Replicon Accession"
chr3_IIIm_head<- merge(chr3_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chr3_IIIm_head<- chr3_IIIm_head%>%arrange(order)


#Chromosome 4
chr4_scaff<- read.csv("chromosome4_scaffold.csv", header = TRUE)
chr4_scaff<- chr4_scaff[c(-1)]
colnames(chr4_scaff)[1]<- "Replicon Accession"
chr4_IIIm_head<- merge(chr4_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chr4_IIIm_head<- chr4_IIIm_head%>%arrange(order)


#Chromosome 5
chr5_scaff<- read.csv("chromosome5_scaffold.csv", header = TRUE)
chr5_scaff<- chr5_scaff[c(-1)]
colnames(chr5_scaff)[1]<- "Replicon Accession"
chr5_IIIm_head<- merge(chr5_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chr5_IIIm_head<- chr5_IIIm_head%>%arrange(order)


#Chromosome X
chrX_scaff<- read.csv("chromosomeX_scaffold.csv", header = TRUE)
chrX_scaff<- chrX_scaff[c(-1)]
colnames(chrX_scaff)[1]<- "Replicon Accession"
chrX_IIIm_head<- merge(chrX_scaff, IIIm_head, by="Replicon Accession")
library(dplyr)
chrX_IIIm_head<- chrX_IIIm_head%>%arrange(order)

write.csv(chr1_IIIm_head, "1stchromosome_hetsIIIm_head.csv")
write.csv(chr2_IIIm_head, "2ndchromosome_hetsIIIm_head.csv")
write.csv(chr3_IIIm_head, "3rdchromosome_hetsIIIm_head.csv")
write.csv(chr4_IIIm_head, "4thchromosome_hetsIIIm_head.csv")
write.csv(chr5_IIIm_head, "5thchromosome_hetsIIIm_head.csv")
write.csv(chrX_IIIm_head, "Xchromosome_hetsIIIm_head.csv")



####### Geetting gene info for the scaffolds. 
# This portion of code should be run on a high performance computing cluster as this takes a long time to do 

## 1st Chromosome

C1<- read.csv("1stchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
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
  for (i in 1:2605){
    for(k in 1:17508){
      if(C1[i,2]==source[k,7]){
        if(C1[i,3]>=source[k,8] & C1[i,3]<=source[k,9]){
          C1[i,9]=source[k,5]
          C1
        }
      }
    }
  }
  return (C1)
}
C1_1<- fun()
write.csv(C1_1,"IIImhets_1stchromosome_with_gene_info.csv")


##2nd chromosome

C2<- read.csv("2ndchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
C2$'Replicon.Accession'<- as.character(C2$'Replicon.Accession')
C2$Replicon.Accession<- gsub("\\..*","",C2$Replicon.Accession)
C2$Position<- as.numeric(C2$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for C2 for gene names
C2$gene = NA
#Getting information from source to file with allele information
fun1<- function(x){
  for (i in 1:4351){
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


### 3rd chromosome

C3<- read.csv("3rdchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
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
  for (i in 1:15517){
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


### 4th chromosome

C4<- read.csv("4thchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
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
  for (i in 1:2364){
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


###5th chromosome

C5<- read.csv("5thchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
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
  for (i in 1:3315){
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


### X chromosome

CX<- read.csv("Xchromosome_hetsIIIm_head.csv", header=TRUE, sep=",")
CX$'Replicon.Accession'<- as.character(CX$'Replicon.Accession')
CX$Replicon.Accession<- gsub("\\..*","",CX$Replicon.Accession)
CX$Position<- as.numeric(CX$Position)
#Reading scaffold and gene information 
source<- read.delim("gene_coord_ME_gb.tsv", header = TRUE)
source$gene<- as.character(source$gene)
#Creating a column for CX for gene names
CX$gene = NA
#Getting information from source to file with allele information
fun5<- function(x){
  for (i in 1:88){
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


###Extracting allele depth from vcf files

ad<- extract.gt(vcf, element = "AD")
ad<- as.data.frame(ad)
ad<- setNames(cbind(rownames(ad), ad, row.names=NULL), c("scaff", "CSrab18Chead", "CSrab29Chead", "IsoCS18Chead", "IsoCS29Chead"))
ad$scaff<- as.character(ad$scaff)
library(dplyr)
ad<- separate(ad, CSrab18Chead, into = c("CSrab18Chead.1", "CSrab18Chead.2"), sep=",")
ad<- separate(ad, CSrab29Chead, into = c("CSrab29Chead.1", "CSrab29Chead.2"), sep=",")
ad<- separate(ad, IsoCS18Chead, into = c("IsoCS18Chead.1", "IsoCS18Chead.2"), sep=",")
ad<- separate(ad, IsoCS29Chead, into = c("IsoCS29Chead.1", "IsoCS29Chead.2"), sep=",")
ad$CSrab18Chead.1<- as.numeric(ad$CSrab18Chead.1)
ad$CSrab18Chead.2<- as.numeric(ad$CSrab18Chead.2)
ad$CSrab29Chead.1<- as.numeric(ad$CSrab29Chead.1)
ad$CSrab29Chead.2<- as.numeric(ad$CSrab29Chead.2)
ad$IsoCS18Chead.1<- as.numeric(ad$IsoCS18Chead.1)
ad$IsoCS18Chead.2<- as.numeric(ad$IsoCS18Chead.2)
ad$IsoCS29Chead.1<- as.numeric(ad$IsoCS29Chead.1)
ad$IsoCS29Chead.2<- as.numeric(ad$IsoCS29Chead.2)

#normalizing allele depths
ad$CSrab18Chead.1<- ad$CSrab18Chead.1*1000000/85784784
ad$CSrab18Chead.2<- ad$CSrab18Chead.2*1000000/85784784
ad$CSrab29Chead.1<- ad$CSrab29Chead.1*1000000/107597223
ad$CSrab29Chead.2<- ad$CSrab29Chead.2*1000000/107597223
ad$IsoCS18Chead.1<- ad$IsoCS18Chead.1*1000000/67102691
ad$IsoCS18Chead.2<- ad$IsoCS18Chead.2*1000000/67102691
ad$IsoCS29Chead.1<- ad$IsoCS29Chead.1*1000000/97019962
ad$IsoCS29Chead.2<- ad$IsoCS29Chead.2*1000000/97019962

#Putting YM together
ad1<- within(ad, IsoCS18Chead<- paste(IsoCS18Chead.1, IsoCS18Chead.2, sep=","))
ad1<- within(ad1, IsoCS29Chead<- paste(IsoCS29Chead.1, IsoCS29Chead.2, sep=","))
ad1<-ad1[, c(1:5,10,11)]
library(reshape2)

ad<- melt(ad, id="scaff")

#For new[1]dataframe
ad1<- melt(ad1, id="scaff")


#### For first chromosome genes

read_IIIm1<- read.csv("IIImhets_1stchromosome_with_gene_info.csv", header=TRUE, sep = ",")

###Reading genes with significant GxT interaction from DESeq2 analysis
interaction<- read.csv("interaction_head.csv", header=TRUE, sep=",")
colnames(interaction)[1]<-"gene"
ASE_IIImhets<- merge(interaction, read_IIIm1, by="gene")  ###No common genes found. So no further analysis was done.


### For second chromosome genes

read_IIIm2<- read.csv("IIImhets_2ndchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIIm2, by="gene")

ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Chead, into = c("CSrab18Chead.1", "CSrab18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Chead, into = c("CSrab29Chead.1", "CSrab29Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Chead, into = c("IsoCS18Chead.1", "IsoCS18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Chead, into = c("IsoCS29Chead.1", "IsoCS29Chead.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Chead<- paste(IsoCS18Chead.1, IsoCS18Chead.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Chead<- paste(IsoCS29Chead.1, IsoCS29Chead.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:16,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:12,]
#ASE_IIImhets1<- ASE_IIImhets1[1:42,]

New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                      "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead",
                                      "IsoCS18Chead.1"="IsoCS18Chead", "IsoCS18Chead.2"="IsoCS18Chead",
                                      "IsoCS29Chead.1"="IsoCS29Chead", "IsoCS29Chead.2"="IsoCS29Chead")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,16,by=8))
    for(k in seq(1,12, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Chead" & New1[k+5,2]=="IsoCS29Chead"){
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
New2<- transform(New2, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                        "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead")))

write.csv(New2, "Heads_2nd_chromosome_heterozygous_data.csv")


### For genes on the third chromosome


read_IIIm3<- read.csv("IIImhets_3rdchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIIm3, by="gene")

ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Chead, into = c("CSrab18Chead.1", "CSrab18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Chead, into = c("CSrab29Chead.1", "CSrab29Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Chead, into = c("IsoCS18Chead.1", "IsoCS18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Chead, into = c("IsoCS29Chead.1", "IsoCS29Chead.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Chead<- paste(IsoCS18Chead.1, IsoCS18Chead.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Chead<- paste(IsoCS29Chead.1, IsoCS29Chead.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:8,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:6,]
#ASE_IIImhets1<- ASE_IIImhets1[1:42,]

New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                      "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead",
                                      "IsoCS18Chead.1"="IsoCS18Chead", "IsoCS18Chead.2"="IsoCS18Chead",
                                      "IsoCS29Chead.1"="IsoCS29Chead", "IsoCS29Chead.2"="IsoCS29Chead")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,8,by=8))
    for(k in seq(1,6, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Chead" & New1[k+5,2]=="IsoCS29Chead"){
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
New2<- transform(New2, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                        "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead")))

write.csv(New2, "Heads_3rd_chromosome_heterozygous_data.csv")


### For fourth chromosome genes

read_IIIm4<- read.csv("IIImhets_4thchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIIm4, by="gene")

ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Chead, into = c("CSrab18Chead.1", "CSrab18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Chead, into = c("CSrab29Chead.1", "CSrab29Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Chead, into = c("IsoCS18Chead.1", "IsoCS18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Chead, into = c("IsoCS29Chead.1", "IsoCS29Chead.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Chead<- paste(IsoCS18Chead.1, IsoCS18Chead.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Chead<- paste(IsoCS29Chead.1, IsoCS29Chead.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:32,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:24,]
#ASE_IIImhets1<- ASE_IIImhets1[1:42,]

New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                      "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead",
                                      "IsoCS18Chead.1"="IsoCS18Chead", "IsoCS18Chead.2"="IsoCS18Chead",
                                      "IsoCS29Chead.1"="IsoCS29Chead", "IsoCS29Chead.2"="IsoCS29Chead")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,32,by=8))
    for(k in seq(1,24, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Chead" & New1[k+5,2]=="IsoCS29Chead"){
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
New2<- transform(New2, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                        "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead")))

write.csv(New2, "Heads_4th_chromosome_heterozygous_data.csv")


### For fifth chromosome genes

read_IIIm5<- read.csv("IIImhets_5thchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIIm5, by="gene")

ASE_IIImhets$scaff<- paste(ASE_IIImhets$Replicon.Accession,".1", sep="")
ASE_IIImhets$scaff<- paste(ASE_IIImhets$scaff,ASE_IIImhets$Position, sep="_")
ASE_IIImhets<- ASE_IIImhets[,c(17,1,12:16)]
library(dplyr)
ASE_IIImhets<- ASE_IIImhets%>%arrange(order)
library(tidyr)
ASE_IIImhets<- separate(ASE_IIImhets, CSrab18Chead, into = c("CSrab18Chead.1", "CSrab18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, CSrab29Chead, into = c("CSrab29Chead.1", "CSrab29Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS18Chead, into = c("IsoCS18Chead.1", "IsoCS18Chead.2"), sep="/")
ASE_IIImhets<- separate(ASE_IIImhets, IsoCS29Chead, into = c("IsoCS29Chead.1", "IsoCS29Chead.2"), sep="/")

#For new[1] dataframe
ASE_IIImhets1<- within(ASE_IIImhets, IsoCS18Chead<- paste(IsoCS18Chead.1, IsoCS18Chead.2, sep=","))
ASE_IIImhets1<- within(ASE_IIImhets1, IsoCS29Chead<- paste(IsoCS29Chead.1, IsoCS29Chead.2, sep=","))
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6, 11:13)]
ASE_IIImhets1<- ASE_IIImhets1[,c(1:6,8:9,7)]

#For first dataframe
library(reshape)
ASE_IIImhets<- melt(ASE_IIImhets, id=c("scaff", "gene"))
ASE_IIImhets<- ASE_IIImhets[1:24,]

#For new[1]dataframe
ASE_IIImhets1<- melt(ASE_IIImhets1, id=c("scaff", "gene"))
ASE_IIImhets1<- ASE_IIImhets1[1:18,]
#ASE_IIImhets1<- ASE_IIImhets1[1:42,]

New<- merge(ASE_IIImhets, ad, by=c("scaff","variable"))

colnames(New)[4]<- "Allele"
colnames(New)[5]<- "Normalized_count"
colnames(New)[2]<- "GT"
New$GT<- as.character(New$GT)
New$scaff<-as.character(New$scaff)
library(plyr)
New<- transform(New, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                      "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead",
                                      "IsoCS18Chead.1"="IsoCS18Chead", "IsoCS18Chead.2"="IsoCS18Chead",
                                      "IsoCS29Chead.1"="IsoCS29Chead", "IsoCS29Chead.2"="IsoCS29Chead")))


#New[1]dataframe
New1<- merge(ASE_IIImhets1, ad1, by=c("scaff", "variable"))
colnames(New1)[4]<- "Allele"
colnames(New1)[5]<- "Normalized_count"
colnames(New1)[2]<- "GT"
New1$GT<- as.character((New1$GT))
New1$scaff<- as.character(New1$scaff)


fun<- function(x){
  for (i in seq(1,24,by=8))
    for(k in seq(1,18, by=6)){
      a=New[i,1]
      if (New[i+4,1]==a & New[i+5,1]==a){
        if (New1[k+4,1]==a & New1[k+5,1]==a){
          if (New[i+6,2]== "IsoCS29Chead" & New1[k+5,2]=="IsoCS29Chead"){
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
New2<- transform(New2, GT=revalue(GT, c("CSrab18Chead.1"="CSrab18Chead","CSrab18Chead.2"="CSrab18Chead",
                                        "CSrab29Chead.1"="CSrab29Chead", "CSrab29Chead.2"="CSrab29Chead")))

write.csv(New2, "Heads_5th_chromosome_heterozygous_data.csv")


### For X chromosome genes


read_IIImX<- read.csv("IIImhets_Xchromosome_with_gene_info.csv", header=TRUE, sep = ",")

ASE_IIImhets<- merge(interaction, read_IIImX, by="gene")   ### No common genes found. So, no further analysis was done



####################### IIIIM-III calculations
#### Chromosome 2


Total<- read.csv("Heads_2nd_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),2)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),2)
Total$order<- 1:12
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:7){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:4){
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
  for (i in 1:1) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:2) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:2){
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

write.csv(df6, "head_2nd_chromosome_genes_details")


#### Chromosome 3


Total<- read.csv("Heads_3rd_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),1)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),1)
Total$order<- 1:6
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:3){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:2){
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

write.csv(df2, "head_3rd_chromosome_genes_details")



### Chromosome 4

Total<- read.csv("Heads_4th_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),4)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),4)
Total$order<- 1:24
Total<- Total[7:12,] #### Keeping only variant sites that have all the necessary information required
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:3){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:2){
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


write.csv(df2, "head_4th_chromosome_genes_details")



### Chromosome 5

Total<- read.csv("Heads_5th_chromosome_heterozygous_data.csv", header=TRUE, sep = ",")
Total$Temp<- rep(c("18C","18C","29C","29C","18C","29C"),3)
Total$allele<- rep(c("III","IIIM","III","IIIM","YM", "YM"),3)
Total$order<- 1:18
df<- Total[!duplicated(Total[c(2,7)]),] #keeping scaffolds and two temperatures
df1<- df[,c(4,2,7)]
df1$'IIIM-III'<- NA
T3<- subset(Total, Total$allele=="III" | Total$allele=="IIIM")


funct<- function(x){
  for (i in 1:11){
    a=T3[i,2]
    b=T3[i,7]
    if (T3[i+1,2]==a & T3[i+1,7]==b){
      for (k in 1:6){
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
  for (i in 1:1) {
    a=df4[i,1]
    df5<- subset(df2, gene==a & Temp=="18C")
    for (k in 1:2) {
      if (df3[k,1]==a & df3[k,2]=="18C"){
        df3[k,3]= mean(df5$IIIM.III)
        df3[k,4]= sd(df5$IIIM.III)/sqrt(length(df5$IIIM.III))
        df3
      }}
    df5<- subset(df2, gene==a & Temp=="29C")
    for (k in 1:2){
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

write.csv(df6, "head_5th_chromosome_genes_details")





