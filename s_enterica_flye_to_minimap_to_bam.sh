#!/bin/bash

#Example flye package installation with bioconda if required. If a non-conda installation is desired, please see https://github.com/mikolmogorov/Flye/blob/flye/docs/INSTALL.md for alternatives.
#conda env create -n flye
#conda activate flye
#conda install flye

#minimap installed via compiling binary executable directly. See https://github.com/lh3/minimap2 for first time installation or alternatives. 

#samtools installed directly with “sudo apt install samtools”. See https://www.htslib.org/download/ for alternatives

#download locations for S.enterica data:
# S.enterica fasta file: https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fasta?acc=SRR32410565
# S.enterica reference genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz

#Actual data download
wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fasta?acc=SRR32410565

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz

gunzip GCF_000006945.2_ASM694v2_genomic.fna.gz


#Genome assembly with flye:
flye --nano-hq SRR32410565.fasta.gz -o Flye_Assembly -t 12

#Genome alignment with minimap. -t12 used due to a 12-core machine for multithreading, substitute as desired.:
minimap2 -ax map-ont -t12 GCA_000006945.2_ASM694v2_genomic.fna Flye_Assembly/assembly.fasta > test1.sam

#File conversion with samtools for
samtools view -bS test1.sam -o test1.bam
samtools sort test1.bam -o test1.sorted.bam
samtools index test1.sorted.bam
