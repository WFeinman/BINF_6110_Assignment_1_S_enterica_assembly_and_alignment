#!/bin/bash

#SRA toolkit required for requested file interpretations, install instructions https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

#Example flye package installation with bioconda if required. If a non-conda installation is desired, please see https://github.com/mikolmogorov/Flye/blob/flye/docs/INSTALL.md for alternatives.
#conda env create -n flye
#conda activate flye
#conda install flye

#minimap installed via compiling binary executable directly. See https://github.com/lh3/minimap2 for first time installation or alternatives. 

#samtools installed directly with “$sudo apt install samtools”. See https://www.htslib.org/download/ for alternatives

#clair3 installed via singularity container: $singularity pull docker://hkubal/clair3:latest. See https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#installation for alternatives

#download locations for S.enterica data:
# S.enterica SRR file:https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565
# S.enterica reference genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz
# S.enterica annotation: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.gtf.gz

#Actual data downloads and unpacking
wget -nc https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

#works, just takes it's sweet time. Like 10-20 minutes.
fasterq-dump SRR32410565

seqtk fqchk SRR32410565.fastq | head

wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz

wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.gtf.gz

gunzip GCA*.gz

#Genome assembly with flye:
flye --nano-hq SRR32410565.fastq -o Flye_Assembly -t 12

seqtk fqchk Flye_Assembly/assembly.fasta | head

#Genome alignment with minimap. -t12 used due to a 12-core machine for multithreading, substitute as desired.:
minimap2 -ax map-ont -t12 \
    GCA_000006945.2_ASM694v2_genomic.fna \
    Flye_Assembly/assembly.fasta > test1.sam

#File conversion with samtools for
samtools view -bS test1.sam -o test1.bam
samtools sort test1.bam -o test1.sorted.bam
samtools index test1.sorted.bam

samtools faidx GCA_000006945.2_ASM694v2_genomic.fna


minimap2 -ax map-ont -t 12 \
  GCA_000006945.2_ASM694v2_genomic.fna \
  SRR32410565.fastq | \
  samtools sort -o test2.sorted.bam && \
  samtools index test2.sorted.bam

#Variant calling with clair3

singularity pull docker://hkubal/clair3:latest

mkdir clair3_out

#note: absolute path required for singularity variant calling, so exporting current pwd for later retrieval
ABSOLUTE_PATH=$(pwd)
export STORED_FILE_PATH=$ABSOLUTE_PATH

singularity exec \
  -B $(echo $STORED_FILE_PATH),$(echo $STORED_FILE_PATH)/clair3_out\
  clair3_latest.sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=$(echo $STORED_FILE_PATH)/test2.sorted.bam \
  --ref_fn=$(echo $STORED_FILE_PATH)/GCA_000006945.2_ASM694v2_genomic.fna \
  --threads=12 \
  --platform="ont" \
  --model_path="/opt/models/r1041_e82_400bps_sup_v500" \
  --include_all_ctgs \
  --output=$(echo $STORED_FILE_PATH)/clair3_out \
  --no_phasing_for_fa \
  --include_all_ctgs \
  --haploid_precise \
  --enable_variant_calling_at_sequence_head_and_tail

gunzip clair3_out/merge_output.vcf.gzm

mv clair3_out/merge_output.vcf merge_output.vcf