## Differential gene expression

### Build hisat2 index

hisat2-build --ss splicesites.txt --exon exons.txt mds_ref_Musca_domestica-2.0.2_names.fa mds_ref_Musca_domestica-Ind 


### Alignment to reference genome


hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Head_S7.fq | samtools view -bS - > CSrab-18C-Head_S7.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Head_S8.fq | samtools view -bS - > CSrab-18C-Head_S8.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Head_S9.fq | samtools view -bS - > CSrab-18C-Head_S9.hisat.bam  
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Head_S1.fq | samtools view -bS - > CSrab-29C-Head_S1.hisat.bam  
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Head_S5.fq | samtools view -bS - > CSrab-29C-Head_S5.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Head_S6.fq | samtools view -bS - > CSrab-29C-Head_S6.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Head_S2.fq | samtools view -bS - > IsoCS-18C-Head_S2.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Head_S3.fq | samtools view -bS - > IsoCS-18C-Head_S3.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Head_S4.fq | samtools view -bS - > IsoCS-18C-Head_S4.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Head_S10.fq | samtools view -bS - > IsoCS-29C-Head_S10.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Head_S11.fq | samtools view -bS - > IsoCS-29C-Head_S11.hisat.bam 
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Head_S12.fq | samtools view -bS - > IsoCS-29C-Head_S12.hisat.bam 


###Sorting bam files


samtools sort -O bam -T tmp_Head1 -o CSrab-29C-H1.sorted.bam CSrab-29C-Head_S1.hisat.bam
samtools sort -O bam -T tmp_Head2 -o CSrab-29C-H2.sorted.bam CSrab-29C-Head_S5.hisat.bam
samtools sort -O bam -T tmp_Head3 -o CSrab-29C-H3.sorted.bam CSrab-29C-Head_S6.hisat.bam

samtools sort -O bam -T tmp_Head4 -o CSrab-18C-H1.sorted.bam CSrab-18C-Head_S7.hisat.bam
samtools sort -O bam -T tmp_Head5 -o CSrab-18C-H2.sorted.bam CSrab-18C-Head_S8.hisat.bam
samtools sort -O bam -T tmp_Head6 -o CSrab-18C-H3.sorted.bam CSrab-18C-Head_S9.hisat.bam

samtools sort -O bam -T tmp_Head7 -o IsoCS-29C-H1.sorted.bam IsoCS-29C-Head_S10.hisat.bam
samtools sort -O bam -T tmp_Head8 -o IsoCS-29C-H2.sorted.bam IsoCS-29C-Head_S11.hisat.bam
samtools sort -O bam -T tmp_Head9 -o IsoCS-29C-H3.sorted.bam IsoCS-29C-Head_S12.hisat.bam

samtools sort -O bam -T tmp_Head10 -o IsoCS-18C-H1.sorted.bam IsoCS-18C-Head_S2.hisat.bam
samtools sort -O bam -T tmp_Head11 -o IsoCS-18C-H2.sorted.bam IsoCS-18C-Head_S3.hisat.bam
samtools sort -O bam -T tmp_Head12 -o IsoCS-18C-H3.sorted.bam IsoCS-18C-Head_S4.hisat.bam



###Counting reads


python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-H1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-H1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-H2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-H2.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-H3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-H3.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-H1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-H1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-H2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-H2.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-H3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-H3.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-H1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-H1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-H2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-H2.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-H3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-H3.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-H1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-H1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-H2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-H2.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-H3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-H3.txt 

