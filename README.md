# NVIDIA-Clara-Parabricks-Somatic-Variant-Calling-on-AWS
## _AWS Blog (Somatic Variant Calling )_

## Data download required to run NVIDIA Clara Parabricks Somatic Variant Calling Pipeline.
### Download publicly available BAM files using wget
### WGS Normal sample (50X)
BAM: wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam
#### File size: 95GB
BAM INDEX: wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bai
### WGS Tumor sample (50X)
BAM: wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bam
#### File size: 107GB
BAM INDEX: wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bai
 
## Download 100x Whole exome dataset from SEQC2 hosted at NCBI:
### WES Normal Sample (80X)
 - wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_N_1.bwa.dedup.bam
 - wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_N_1.bwa.dedup.bai

### WES Tumor Sample (100X)
 - wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_T_1.bwa.dedup.bam
 - wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_T_1.bwa.dedup.bai
 
## Downloading GRCh38 reference genome along with known sites from GDC/SEQC2, these files will be used for alignment and variant calling.
The reference files are hosted at GDC/NCI website. (https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
###	Fasta reference files: GRCh38.d1.vd1.fa.tar.gz
    wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz
    tar -xvzf GRCh38.d1.vd1.fa.tar.gz
###	BWA index reference files: GRCh38.d1.vd1_BWA.tar.gz
    wget https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225 -O GRCh38.d1.vd1_BWA.tar.gz
    tar -xvzf GRCh38.d1.vd1_BWA.tar.gz 
    
### dbSNP file to use : dbSNP146
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
### InDel file from 1k genome to use for known sites : 1000 genome Phase 3
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
    ## Please edit line 22 of this file add missing double quote before dbSNP
       Description=“dbSNP ssID of the allele” 
       ##INFO=<ID=ssID,Number=A,Type=String,Description=dbSNP ssID of the allele">
    ## Please edit line 42 change “POS=POS-1” to “POS_POS-1”
       ##INFO=<ID=POS=POS1,Number=0,Type=Flag,Description="POS has been adjusted due to missing REF in NCBI VCF file">
### Curated InDel file from: Mills_and_1000G_gold_standard.indels
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz

##	Before starting the analysis run, please install following software’s will be handy for some data processing:Download accessory tools 
```
  cd /mnt/disks/local
  ## Downlaod and install GATK 
  wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
  unzip gatk-4.2.0.0.zip
  ## Install samtools
  sudo apt install samtools
  ## Install Tabix
  sudo apt install tabix
  ## index fasta file
  samtools faidx Refs/GRCh38.d1.vd1.fa
  ## create dictionary file from fasta: 
  gatk-4.2.0.0/gatk CreateSequenceDictionary R=/home/ubuntu/Refs/GRCh38.d1.vd1.fa O=/home/ubuntu/Refs/GRCh38.d1.vd1.dict
```
