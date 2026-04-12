# BINF_6110_Assignment_1_S_enterica_assembly_and_alignment

## Overal Goal:
The script aims to take a case study in the form of FastQ reads from an Oxford Nanopore R10 sequence for a Salmonella enterica sample. These reads are to be assembled into a consensus genome, compared to a reference genome, and have any variants visualized as both an assembly and a functional; annotation alignment. The code to do so is assembled in a functional demonstration vignette as an executable shell script in this repository.


## Shell Script Execution:
Ensure flye, FastQC, minimap2, and samtools are properly installed and active. Tested on Windows Subsystem for Linux version 2.5.10.0. Last confirmed functional package versions: flye v2.9.6-b1802, FastQC v0.12.1, minimap2 v2.28-r1209, samtools v1.23.1. 
  - Place the “s_enterica_flye_to_minimap_to_bam.sh” script in desired file location. 
- Give the script execution privileges with:
-     chmod +x s_enterica_flye_to_minimap_to_bam.sh

  And for repeated scripts:
-	  chmod +x s_enterica_no_fastq_download_flye_to_minimap_to_bam.sh

Execute the priary script with:
-     bash s_enterica_flye_to_minimap_to_bam.sh

On execution, requisite genomic data and reference genome will be downloaded to the working directory.  

Alternatively, should the file `SRR32410565.fastq` already be downloaded from a previous script, an optional script without that specific download has been included for ease of repetition: 
-     bash ss_enterica_no_fastq_download_flye_to_minimap_to_bam.sh

Either way, the above genome assembly and alignment will be conducted on the genomic data using the downloaded reference. 

Once completed, a `test1.sorted.bam` file will be generated in the working directory. This may be loaded into the IGV and compared to the reference genome `GCF_000006945.2_ASM694v2_genomic.fna` for the assembly. 


## Visualization:
This data is intended to be viewed with the Integrative Genomics Viewer (IGV) for Desktop, downloaded seperately here: `https://igv.org/doc/desktop/#DownloadPage/` Last tested for IGV Windows Desktop version 2.19.7.

In addition to the assemblyfile `test1.sorted.bam`, the `test2.sorted.bam` file may be likewise compared to the reference for the aligned sequences.

The `GCA_000006945.2_ASM694v2_genomic.gtf` file can be loaded as a broad overview of known gene functions within S. enterica.

Finally, the `merge_output.vcf` highlights shared SNPs, indels, and other variants where the sample reads as a whole significant;y differ from the reference genome.

With the Integrative Genomics Viewer, all of these features may be compared to the ference genome seperately, or at the same time.

