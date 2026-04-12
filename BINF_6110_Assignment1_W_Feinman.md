# BINF 6110 Assignment 1:

## Demonstration of *S. Enterica* Genome Assembly, Annotation, and Alignment.

William Feinman

**[Background:]{.underline}**

Genome assembly is a vital process to modern genetic analysis pipelines.
For most forms of genetic sequencing, a large number of DNA reads
(nucleotide sequences) of unclear origin and variable length are
produced. In order to make use of this data, genome assembly is required
to sort these reads, determine their origin, and map these reads in
relation to each other. This is especially important for use in aligning
reads to other genomes, where proper sequence and orientation can be
vital for determining identity and finding differences. Depending on the
species, sample source, sequencing technology, and intended use, the
best genome assembly tool for a given analysis can vary significantly.
When the subject of concern is pathogenic, this is emphasized all the
further. (Bogaerts et al., 2025)

The main differences in genetic analysis approaches begin with the type
of sequencing technology used to produce the reads. With the exception
of short RNA sequences and viroids, genomes and DNA/RNA strands were
typically too long to accurately analyze as a single molecule in a
practical timeframe until the refinement of long-read sequencing. As
such, whole genome sequencing had been performed using fragment based
short-read approaches like Illumina sequencing. In such sequencing, the
tested genome is intentionally cut into short sequences of a couple of
hundred nucleotides in length (typically, 150-200nt, though wider ranges
exist) to benefit from parallel sequencing. Overlapping motifs between
fragments are used to reconstruct a model of the original genome via
assembly.

This short read sequencing, while accurate, has a problem assembling the
fragments into a larger architecture, or dealing with repeated
sequences. Contiguous reads (contigs), while allowing fragments to be
put together in larger groups, may not be able to overlap into a single
unified whole, forcing assembly software to guess at cross-contig order.
If motifs within a read are repeated for longer than the length of that
read, it may not be possible to accurately tell how many of the reads
are due to duplicate fragments as an artifact of sequencing and how many
are due to genuine repeats in the original genome.

This challenge has knock-on effects for later alignment to a reference
genome, which is vital for comparison between samples. To be useful for
comparative analysis, a sample sequence must be aligned against a known
sample, and differences compared. Without an accurate repeat count or
genomic architecture, major sources of variation between individuals and
species can be overlooked, misattributed, or mis-identified. For
example, a translocation event may be missed entirely if the relevant
gene was only found in an isolated contig without sufficient surrounding
context sequenced.

Long read sequencing, on the other hand, attempts to preserve DNA/RNA
strands intact and sequence them as one unit, using technological
advances such as electric-resistance-reporting synthetic membrane-bound
proteins (Oxford Nanopore) to sequence nucleotide strands far faster
than traditional methods. As a consequence of the long read length,
overlap sequences between strands are far more likely, allowing the
sequence to be assembled into more complete, less fragmented contigs.

Earlier in their development, however, this increased architectural
consistency came at a cost in individual nucleotide accuracy. The longer
a single read, the higher the chance of accumulating error from a
variety of factors (In recent years, however, long read sequencing such
as Oxford Nanopore has grown sufficiently accurate to be comparable to
short read sequencing while addressing these architecture issues.

This paper aims to outline common goals and challenges of assembling a
genome and aligning it to a reference genome. Then, using that
information, the script aims to take a case study in the form of FastQ
reads from an Oxford Nanopore R10 sequence for a *Salmonella enterica*
sample (Bogaerts et al., 2025), as an example of both interesting
genomics and pathogenic importance.

These reads are to be assembled into a consensus genome, compared to a
reference genome, and have any variants visualized as both an assembly
and a functional; annotation alignment with multiple useful metrics of
comparison to highlight variants of interest in the sequenced sample.
The code to do so is assembled in a functional demonstration vignette as
an executable bash script in the attached github repository, replicated
below for sake of completeness:
<https://github.com/WFeinman/BINF_6110_Assignment_1_S_enterica_assembly_and_alignment>

**[Methods:]{.underline}**

Sequence quality was confirmed by running FastQC on the input fastq
file. Though specialized for short reads, FastQC is a way of flagging
pressing issues in a dataset, such as a low quality score or bad reads.
Here, it identified done. Genome assembly was performed using flye
version 2.9.6-b1802. Flye is a longstanding, well documented genome
annotation software with simple installation and quick operation, making
it highly beneficial as a learning tool. Additionally, it is specialized
as a long-read assembly tool with explicit argument support for R10
Oxford Nanopore reads such as ours. Other assembly software may be more
accurate, but take longer and are more technically involved. Unicycle in
particular, though a recent paper showed it had very promising
performance metrics, has two other competing assembly software programs
put out by the same developer, and is in maintenance rather than active
development.

Genome alignment was performed using minimap2 version 2.28-r1209. File
conversion with samtools version 1.23.1. Subsequent genome annotation
was performed with Clair3 v2.0.0 on the un-assembled alignment, as it
has beenone of the top performing annotation packages in recent years
(Helal et al., 2022). Finally .bam files were viewed with Integrative
Genomics Viewer for Windows Desktop version 2.19.7, comparing to
original reference genome from GenBank for S. Enterica.

S. enterica genomic data and Genbank reference genome download locations
included in script file and github readme. All operations performed on
Windows Subsystem for Linux version 2.5.10.0 bash shell, save for final
viewing on the Integrative Genomics Viewer (IGV) for Windows version
2.19.7.

Notably, the long-read dataset has been pre-selected for high read
quality and no exceptionally long reads (expected accuracy score Q20+,
N50:5-15kb), mitigating potential assembly inaccuracies and falling well
within the working range of our alignment software.

If working with lower expected accuracy such as from an earlier
Oxford-nanopore set, the assembly portion of the script should be
adjusted accordingly, such as by using the "--nano-raw" argument for the
flye command in place of the "--nano-hq" to indicate a lower expected
accuracy. If there were a few reads in particular of lower quality than
the rest, they could be filtered out with a quality score rule to
preserve the rest of the data.

**[Bash Script Execution:]{.underline}**

Ensure flye, FastQC, minimap2, and samtools are properly installed and
active. Place the "s_enterica_flye_to_minimap_to_bam.sh" script in
desired file location.

Give the script execution privileges with: chmod +x
s_enterica_flye_to_minimap_to_bam.sh

-   And for repeated scripts:

chmod +x s_enterica_no_fastq_download_flye_to_minimap_to_bam.sh

Execute the script with: bash s_enterica_flye_to_minimap_to_bam.sh

On execution, requisite genomic data and reference genome will be
downloaded to the working directory. Alternatively, should the file
SRR32410565.fastq already be downloaded from a previous script, an
optional script without that specific download has been included for
ease of repetition:

bash s_enterica_no_fastq_download_flye_to_minimap_to_bam.sh

Either way, the above genome assembly and alignment will be conducted on
the genomic data using the downloaded reference.

Once completed, a "test1.sorted.bam" file will be generated in the
working directory. This may be loaded into the IGV and compared to the
reference genome "GCF_000006945.2_ASM694v2_genomic.fna" for the
assembly. The "test2.sorted.bam" file may be likewise compared to the
reference for the aligned annotation.

**[Primary Script file contents:]{.underline}**

#!/bin/bash

#SRA toolkit required for requested file interpretations, install
instructions
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

#Example flye package installation with bioconda if required. If a
non-conda installation is desired, please see
https://github.com/mikolmogorov/Flye/blob/flye/docs/INSTALL.md for
alternatives.

#conda env create -n flye

#conda activate flye

#conda install flye

#minimap installed via compiling binary executable directly. See
https://github.com/lh3/minimap2 for first time installation or
alternatives.

#samtools installed directly with "\$sudo apt install samtools". See
https://www.htslib.org/download/ for alternatives

#clair3 installed via singularity container: \$singularity pull
docker://hkubal/clair3:latest. See
https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#installation for
alternatives

#download locations for S.enterica data:

\# S.enterica SRR
file:https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

\# S.enterica reference genome:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz

\# S.enterica annotation:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.gtf.gz

#Actual data downloads and unpacking

wget -nc
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565

#works, just takes it\'s sweet time. Like 10-20 minutes.

fasterq-dump SRR32410565

seqtk fqchk SRR32410565.fastq \| head

wget -nc
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz

wget -nc
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.gtf.gz

gunzip GCA\*.gz

#Genome assembly with flye:

flye \--nano-hq SRR32410565.fastq -o Flye_Assembly -t 12

seqtk fqchk Flye_Assembly/assembly.fasta \| head

#Genome alignment with minimap. -t12 used due to a 12-core machine for
multithreading, substitute as desired.:

minimap2 -ax map-ont -t12 \\

GCA_000006945.2_ASM694v2_genomic.fna \\

Flye_Assembly/assembly.fasta \> test1.sam

#File conversion with samtools for

samtools view -bS test1.sam -o test1.bam

samtools sort test1.bam -o test1.sorted.bam

samtools index test1.sorted.bam

samtools faidx GCA_000006945.2_ASM694v2_genomic.fna

minimap2 -ax map-ont -t 12 \\

GCA_000006945.2_ASM694v2_genomic.fna \\

SRR32410565.fastq \| \\

samtools sort -o test2.sorted.bam && \\

samtools index test2.sorted.bam

#Variant calling with clair3

singularity pull docker://hkubal/clair3:latest

mkdir clair3_out

#note: absolute path required for singularity variant calling, so
exporting current pwd for later retrieval

ABSOLUTE_PATH=\$(pwd)

export STORED_FILE_PATH=\$ABSOLUTE_PATH

singularity exec \\

-B \$(echo \$STORED_FILE_PATH),\$(echo \$STORED_FILE_PATH)/clair3_out\\

clair3_latest.sif \\

/opt/bin/run_clair3.sh \\

\--bam_fn=\$(echo \$STORED_FILE_PATH)/test2.sorted.bam \\

\--ref_fn=\$(echo
\$STORED_FILE_PATH)/GCA_000006945.2_ASM694v2_genomic.fna \\

\--threads=12 \\

\--platform=\"ont\" \\

\--model_path=\"/opt/models/r1041_e82_400bps_sup_v500\" \\

\--include_all_ctgs \\

\--output=\$(echo \$STORED_FILE_PATH)/clair3_out \\

\--no_phasing_for_fa \\

\--include_all_ctgs \\

\--haploid_precise \\

\--enable_variant_calling_at_sequence_head_and_tail

gunzip clair3_out/merge_output.vcf.gzm

mv clair3_out/merge_output.vcf merge_output.vcf

**[Results:]{.underline}**

![](media/image1.png){width="3.361920384951881in"
height="3.3040912073490816in"}

**Figure 1:** Read quality was confirmed with FastQC. Phred score well
above 20 indicates accurate reads.

![](media/image2.png){width="5.984476159230097in"
height="7.50680883639545in"}

**Figure 2:** Example cross-file comparison chart. Top-mid histogram
thickness indicates extent of coverage in region. Color change on same
indicates variance from reference. Red, green, blue, or orange indicates
T, A, C, G substitution (respectively) in place of the reference
nucleotide. Purple indicates insertion. White indicates a read gap.
Example of a low coverage region being highly variable from reference.

![](media/image3.png){width="1.740961286089239in"
height="6.096423884514436in"}
![](media/image4.png){width="1.7028357392825897in"
height="6.2429319772528435in"}
![](media/image5.png){width="3.0428346456692914in"
height="6.165593832020997in"}

**Figure 3:** Example of distinct chimeric alignment. Read does not map
cleanly to aligned region, so supplementary alignment recorded at the
closest match site.In this case, a high coverage regtion shows high
variation from the reference.

![](media/image6.png){width="6.5in" height="6.370138888888889in"}

**Figure 4:** Example of secondary alignment. Unlike chimeric
alignments, reads clearly map to aligned reference region, but a
secondary motif is actually more common than the reference sequnce at a
given spot. In this case, a conserved insertion sequence

![](media/image7.png){width="6.5in" height="5.102777777777778in"}

**Figure 5:** Example of a Salmonella pathogenecity island (SPI-1).
Despite an absence of clear variant calls like prior examples, many
reads show increased insertions in this broad area.

**[Discussion:]{.underline}**

Overall, the *S. enterica* daataset showed a remarkable number of
illustrative alignment features, serving as a good illustrative example
of some of the challenges in categorizing genomes even with and accurate
dataset. (SAM 2025). Even with high quality long reads, low coverage of
a region can lead to interpretation difficulty. Spotty coverage may lead
to important regions being under-identified. Disruption from chimeric or
secondary alignments can allow important traits to be mis-identified.
*Salmonella* species are a highly illustrative example of this. (Sia et
al., 2025)

Such species are economically devastating, highly destructive to animal
species, and have proven resilient to medical intervention (Jajere,
2019). Plasmid swapping, co-infection, and mutable genomes provide
substantially increase variability and evolutionary flexibility, novel
functions being gained and lost through both conventional and
unconventional methods of gene transfer. (Tautz & Domazet-Lošo, 2011)

*Salmonella*'s pathogenicity islands (SPIs) are an illustration of this
behaviour in action. These are genetic regions shared across
*Salmonella* species encoding traits which increase bacterial virulence
and harm to their host. SPI-1, for example, designates a collection of
genes which code for an injectisome protein complex (SipB, SipC, SipD)
aiding in host cell invasion. Other SPIs might affix more firmly in host
guts, evade immune response, increase toxicity, or a variety of other
deleterious effects. (Lou et al., 2019) (Wood et al., 2025)

These regions are not always, however, easily identifiable and
trackable. As with other parts of the *Salmonella* genome, these motifs
may duplicate, split, translocate, and otherwise be modified and move
wildly across genomes and species. These genetic effects provide
challenges to traditional bioinformatic analysis, stressing the
importance of improving existing tools and developing new ones. Assembly
may struggle when gene location and content is highly variable even
within the same species. Alignment is difficult in such cases even with
deep, accurate sequencing, as the genetic relations they are trying to
capture are often non linear, such as horizontal gene transfer.

It should be emphasize that what is captured here was a particularly
deep, accurate sequencing dataset and extensive coverage, which lets us
characterize the gene interactions here with much more confidence. In
studies lacking the time or funding necessary to obtain that accuracy
and depth, capturing these shifts becomes far more difficult. The
importance of a robust genetic analysis pipeline is emphasized,
expecially in light of the last pandemic.

**[\
]{.underline}**

**[References:]{.underline}**

Bogaerts, B., Maex, M., Commans, F., Goeders, N., Van Den Bossche, A.,
De Keersmaecker, S. C. J., Roosens, N. H. C., Ceyssens, P.-J., Mattheus,
W., & Vanneste, K. (2025). Oxford Nanopore Technologies R10 sequencing
enables accurate cgMLST-based bacterial outbreak investigation of
Neisseria meningitidis and Salmonella enterica when accounting for
methylation-related errors. Journal of Clinical Microbiology, 63(10),
e00410-25. <https://doi.org/10.1128/jcm.00410-25>

-   Direct Download:
    <https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565>

-   Reference Genome:
    <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz>

Helal, A. A., Saad, B. T., Saad, M. T., Mosaad, G. S., & Aboshanab, K.
M. (2022). Evaluation of the Available Variant Calling Tools for Oxford
Nanopore Sequencing in Breast Cancer. Genes, 13(9), 1583.
<https://doi.org/10.3390/genes13091583>

Jajere, S. M. (2019). A review of Salmonella enterica with particular
focus on the pathogenicity and virulence factors, host specificity and
antimicrobial resistance including multidrug resistance. Veterinary
World, 12(4), 504--521. <https://doi.org/10.14202/vetworld.2019.504-521>

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of
long, error-prone reads using repeat graphs. Nature Biotechnology,
37(5), 540--546. https://doi.org/10.1038/s41587-019-0072-8

Lou, L., Zhang, P., Piao, R., & Wang, Y. (2019). Salmonella
Pathogenicity Island 1 (SPI-1) and Its Complex Regulatory Network.
Frontiers in Cellular and Infection Microbiology, 9, 270.
<https://doi.org/10.3389/fcimb.2019.00270>

NCBI SRA. (2026). National Center for Biotechnology Information.
<https://github.com/ncbi/sra-tools>

SAM/BAM Format Specification Working Group (2025). Sequence
Alignment/Map Format Specification. Global Alliance for Genomics &
Health. <https://samtools.github.io/hts-specs/SAMv1.pdf>

Sia, C. M., Pearson, J. S., Howden, B. P., Williamson, D. A., & Ingle,
D. J. (2025). Salmonella pathogenicity islands in the genomic era.
Trends in Microbiology, 33(7), 752--764.
<https://doi.org/10.1016/j.tim.2025.02.007>

Tautz, D., & Domazet-Lošo, T. (2011). The evolutionary origin of orphan
genes. Nature Reviews Genetics, 12(10), 692--702.
<https://doi.org/10.1038/nrg3053>

Wood, G., Johnson, R., Powell, J., Bryant, O. J., Lastovka, F., Brember,
M. P., Tourlomousis, P., Carr, J. P., Bryant, C. E., & Chung, B. Y. W.
(2025). The Salmonella pathogenicity island 1 injectisome reprograms
host cell translation to evade the inflammatory response. Nature
Communications, 16(1), 9742.
<https://doi.org/10.1038/s41467-025-64744-w>
