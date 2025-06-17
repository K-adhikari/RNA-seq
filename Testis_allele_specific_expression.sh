##Use STAR to align reads to reference genome

## Generate reference genome files for alignment

STAR --runMode genomeGenerate --genomeDir  --genomeFastaFiles mds_ref_Musca_domestica-2.0.2_names.fa --runThreadN 4 --limitGenomeGenerateRAM 64000000000 > STARgenomeGenerate.out 


## Perform first pass of alignment

mkdir 1pass
cd 1pass
STAR --genomeDir --readFilesIn CSrab-18C-Testes_S2.fq,CSrab-18C-Testes_S3.fq,CSrab-18C-Testes_S4.fq --runThreadN 4 > 1pass.out


## Create a new index using the splice junction information from the first pass

mkdir 2pass
cd 2pass
STAR --runMode genomeGenerate --genomeDir /2pass/ --genomeFastaFiles mds_ref_Musca_domestica-2.0.2_names.fa --sjdbFileChrStartEnd 1pass/SJ.out.tab --sjdbOverhang 99 --genomeChrBinNbits 18 --runThreadN 4 > ./2pass.out


## Perform final pass of alignment

mkdir Finalpass
cd Finalpass
STAR --genomeDir /2pass/ --readFilesIn CSrab-18C-Testes_S2.fq,CSrab-18C-Testes_S3.fq,CSrab-18C-Testes_S4.fq --runThreadN 4 > Final Finalpass.out


## Add read group information, sort, mark duplicates and index

java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ picard-tools-1.133/picard.jar AddOrReplaceReadGroups I=FinalPass/Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=1 RGLB=library RGPL=illumina RGPU=machine RGSM=CSrab18Ctestis TMP_DIR=/tmp/ > AddOrReplaceReadGroups.out 

java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ picard-tools-1.133/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics TMP_DIR=/tmp > MarkDuplicates.out


## SPlit'N'Trim and reassign mapping qualities

java -Xmx8g -jar -Djava.io.tmpdir= /tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T SplitNCigarReads -R mds_ref_Musca_domestica-2.0.2_names.fa -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS > SplitNTrim.out


## Indel realignment

java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R mds_ref_Musca_domestica-2.0.2_names.fa -I split.bam -o target_intervals.list > RealignerTargetCreator.out

java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T IndelRealigner -R mds_ref_Musca_domestica-2.0.2_names.fa -I split.bam -targetIntervals target_intervals.list -o realigned_reads.bam > IndelRealigner.out 


## Filter highest confidence SNPs and INDELs
## We used variant calls from previous RNA-seq analysis to recalibrate reads

#1. Analyze patterns of covariation in the sequence dataset
java -Xmx8g -jar -Djava.io.tmpdir= /tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T BaseRecalibrator -R mds_ref_Musca_domestica-2.0.2_names.fa -I realigned_reads.bam -knownSites CS-B10_JointGenotype_RefSeq_sorted.vcf -o recal_data.table > recal_data.out

#2. Do a second pass to analyze covariation remaining after recalibration
java -Xmx8g -jar -Djava.io.tmpdir= /tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T BaseRecalibrator -R mds_ref_Musca_domestica-2.0.2_names.fa -I realigned_reads.bam -knownSites CS-B10_JointGenotype_RefSeq_sorted.vcf -BQSR recal_data.table -o post_recal.table > post_recal.out 

#3. Generate before/after plots
java -Xmx8g -jar -Djava.io.tmpdir= /tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T AnalyzeCovariates -R mds_ref_Musca_domestica-2.0.2_names.fa -before recal_data.table -after post_recal.table -plots recal_plots.pdf > recal_plots.out

#4. Apply the recalibration to sequence data
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T PrintReads -R mds_ref_Musca_domestica-2.0.2_names.fa -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam > recal_reads.out


## Final variant calling and filtering on recal_reads.bam
# Variant calling
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R mds_ref_Musca_domestica-2.0.2_names.fa -I recal_reads.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o CSrab18Ctestis_VariantCall.vcf > CSrab18Ctestis_VariantCall.out

# Variant filtering
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T VariantFiltration -R mds_ref_Musca_domestica-2.0.2_names.fa -V CSrab18Ctestis_VariantCall.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o CSrab18Ctestis_VariantFilter.vcf > CSrab18Ctestis_VariantFilter.out 

# Output records for joint genotyping
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R mds_ref_Musca_domestica-2.0.2_names.fa -I recal_reads.bam -o CSrab18Ctestis_VariantCall.g.vcf -ERC GVCF > CSrab29Ctestis_VariantCall_gvcf.out 

######################
### We repeated all of the above steps for CSrab at 29°C, IsoCS at 18°C, and IsoCS at 29°C and we used the g.vcf files for joint genotyping
######################


## Genotype
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R mds_ref_Musca_domestica-2.0.2_names.fa -V CSrab18Ctestis_VariantCall.g.vcf -V CSrab29Ctestis_VariantCall.g.vcf -V IsoCS18Ctestis_VariantCall.g.vcf -V IsoCS29Ctestis_VariantCall.g.vcf -stand_call_conf 20.0 -stand_emit_conf 20.0 -o Testis_JointGenotype.vcf > Testis_JointGenotype.out

## Filter variants from Joint Genotyping 
java -Xmx8g -jar -Djava.io.tmpdir=/tmp/ GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T VariantFiltration -R mds_ref_Musca_domestica-2.0.2_names.fa -V Testis_JointGenotype.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Testis_JointGenotype_Filtered.vcf > Testis_JointGenotype_Filtered.out

##Writing filtered variants into a new csv file
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' Testis_JointGenotype_Filtered.vcf > Testis_JointGenotype_Filtered_Pass.vcf