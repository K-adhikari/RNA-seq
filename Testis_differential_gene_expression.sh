### Alignment to reference genome
### Use the same hisat2 index built for head data analysis

hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Testes_S2.fq | samtools view -bS - > CSrab-18C-Testes_S2.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Testes_S3.fq | samtools view -bS - > CSrab-18C-Testes_S3.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-18C-Testes_S4.fq | samtools view -bS - > CSrab-18C-Testes_S4.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Testes_S1.fq | samtools view -bS - > CSrab-29C-Testes_S1.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Testes_S5.fq | samtools view -bS - > CSrab-29C-Testes_S5.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U CSrab-29C-Testes_S6.fq | samtools view -bS - > CSrab-29C-Testes_S6.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Testes_S10.fq | samtools view -bS - > IsoCS-18C-Testes_S10.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Testes_S11.fq | samtools view -bS - > IsoCS-18C-Testes_S11.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-18C-Testes_S12.fq | samtools view -bS - > IsoCS-18C-Testes_S12.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Testes_S7.fq | samtools view -bS - > IsoCS-29C-Testes_S7.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Testes_S8.fq | samtools view -bS - > IsoCS-29C-Testes_S8.hisat.bam
hisat2 -x mds_ref_Musca_domestica-Ind -p 8 --dta-cufflinks -U IsoCS-29C-Testes_S9.fq | samtools view -bS - > IsoCS-29C-Testes_S9.hisat.bam


### Sorting bam files

samtools sort -O bam -T tmp_Testes1 -o CSrab-29C-T1.sorted.bam CSrab-29C-Testes_S1.hisat.bam
samtools sort -O bam -T tmp_Testes2 -o CSrab-29C-T2.sorted.bam CSrab-29C-Testes_S5.hisat.bam
samtools sort -O bam -T tmp_Testes3 -o CSrab-29C-T3.sorted.bam CSrab-29C-Testes_S6.hisat.bam

samtools sort -O bam -T tmp_Testes4 -o CSrab-18C-T1.sorted.bam CSrab-18C-Testes_S2.hisat.bam
samtools sort -O bam -T tmp_Testes5 -o CSrab-18C-T2.sorted.bam CSrab-18C-Testes_S3.hisat.bam
samtools sort -O bam -T tmp_Testes6 -o CSrab-18C-T3.sorted.bam CSrab-18C-Testes_S4.hisat.bam

samtools sort -O bam -T tmp_Testes7 -o IsoCS-29C-T1.sorted.bam IsoCS-29C-Testes_S7.hisat.bam
samtools sort -O bam -T tmp_Testes8 -o IsoCS-29C-T2.sorted.bam IsoCS-29C-Testes_S8.hisat.bam
samtools sort -O bam -T tmp_Testes9 -o IsoCS-29C-T3.sorted.bam IsoCS-29C-Testes_S9.hisat.bam

samtools sort -O bam -T tmp_Testes10 -o IsoCS-18C-T1.sorted.bam IsoCS-18C-Testes_S10.hisat.bam
samtools sort -O bam -T tmp_Testes11 -o IsoCS-18C-T2.sorted.bam IsoCS-18C-Testes_S11.hisat.bam
samtools sort -O bam -T tmp_Testes12 -o IsoCS-18C-T3.sorted.bam IsoCS-18C-Testes_S12.hisat.bam


### Counting reads

python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-T1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-T1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-T2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-T2.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-18C-T3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-18C-T3.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-T1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-T1.txt 
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-T2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-T2.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id CSrab-29C-T3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_CSrab-29C-T3.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-T1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-T1.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-T2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-T2.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-18C-T3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-18C-T3.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-T1.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-T1.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-T2.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-T2.txt
python -m HTSeq.scripts.count -f bam -s reverse -i gene_id IsoCS-29C-T3.sorted.bam ref_Musca_domestica-2.0.2_top_level.gtf > htseq_IsoCS-29C-T3.txt