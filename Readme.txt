This repository consists of codes for analyzing RNA-seq data generated from head and testis tissues of housefly. The raw RNA-seq data are available from NCBI Gene Expression Omnibus accession GSE136188. Additional files used for analysis can be accessed from Texas Data Repository under https://doi.org/10.18738/T8/D14TGI. Associated paper can be accessed from https://onlinelibrary.wiley.com/doi/10.1111/mec.16148

For gene expression study, the reads were mapped to the reference genome using HISAT2 with default settings. The aligned reads were sorted using SAMtools and the reads aligning to genes were counted using HTSeq count. Only the uniquely mapped reads were used, whereas the reads with ambiguous mapping and reads with mapping quality of less than 10 were excluded from the analysis. DESeq2 package in R was used for differential gene expression.

For allele specific expression study was carried out by following GATK best practices workflow for single nucleotide polymorphism (SNP) and indel calling.

