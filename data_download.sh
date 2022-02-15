#!/bin/bash

# Download Normal Bam
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam  
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bai

# Download Tumor Bam
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bai  


# Download Reference
wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz
tar -xvzf GRCh38.d1.vd1.fa.tar.gz
wget https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225 -O Refs/GRCh38.d1.vd1_BWA.tar.gz
tar -xvzf GRCh38.d1.vd1_BWA.tar.gz 


wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
# please check the Readme to edit the below file once downloaded
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz


wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
unzip gatk-4.2.0.0.zip

sudo apt install samtools
sudo apt install tabix
samtools faidx GRCh38.d1.vd1.fa
gatk-4.2.0.0/gatk CreateSequenceDictionary R=GRCh38.d1.vd1.fa O=GRCh38.d1.vd1.dict
