---
Link to files: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3256650.svg)](https://doi.org/10.5281/zenodo.3256650)


questions:
    - TODO
    - What are histone modifications?
    - How is a raw set of scChIC-seq data processed?
    - How do levels in H3K4me3 histone modification differ across cell types in the bone marrow?
objectives:
    - TODO
    - Inspect the read quality
    - Trim low quality bases
    - Map reads on a reference genome
    - Assess the quality of a scChIC-seq experiment
    - Extract coverage files
    - Call enriched regions or peaks
time_estimation: "3h"
key_points:
    - TODO: Key point 1
    - Key point 2
    - Key point 3
contributors:
    - Galaxy Training Team
    - vivekbhr
    - jakeyeung
    - anna-alemany
    - BuysDB
    - Alexander van Oudenaarden
software needed:
    - fastqc
    - samtools
    - deeptools
    - bwa
    - hiddenDomains
    - macs2
    - R
    - perl
    - python




# Introduction

Within a cell nucleus, the DNA is tightly-packed and the chromatin is spatially distributed with different levels and scales of organizations. At the smallest scale, DNA is packaged into units called nucleosomes, made of eight histone proteins.


Transcription factors (TFs) in concert with histone modifications shape the chromatin landscape of the genome, and thus regulate cell types and cell states. Histone modifications form an adaptable epigenetic regulatory layer that mediate dynamic transcriptional programs. Functional genomics assays, the most popular involving chromatin immunoprecipitation (ChIP), have revealed active and repressive chromatin structures in bulk tissues. However, inefficiencies of ChIP hinder its application in single cells, preventing genome-wide analysis of histone modifications along the continuum of cellular phenotypes. Therefore, how chromatin landscapes change between repressed, poised, and active states during development and homeostasis is relatively unexplored at the single-cell level.

Binding certain proteins to each of the eight histone proteins may modify the chromatin structure and may result in changes in transcription level. For example, the H3K4me3 is adding 3 methyl-group of the 4th Lysine in the histone 3 amino-acid. This modification is known to activate the transcription on nearby genes by opening the chromatin. The H3K27me3 on the other hand is inactivating the transcription of the nearby genes:

![Fadloun et al, 2013](images/histone_modifications.jpg "Source: Fadloun et al, 2013")


In the upcoming tutorial, we will look at H3K4me3 scChIC-seq data from mouse bone marrow and see how we can do preliminary analyses to see how histone modification levels differ between cell types. We have prepared the files such that the individual cells are grouped into three clusters based on three cell types: erythroblasts, lymphocytes, and granulocytes. Your job is to infer which cluster corresponds to which cell type.


- TODO: Add Sample names here


# Step 1: Quality control and treatment of the sequences

Let's first assess the quality of the demultiplexed `fastq` files. Demultiplexed means UMI and cell barcodes have been clipped from the `fastq` file and placed in the header.

> ### Hands-on: First look at the `fastq` files
> 1. Import `demultiplexedR1_10000rows.fastq.gz` and `demultiplexedR2_10000rows.fastq.gz` from [Zenodo](https://zenodo.org/record/1324070) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/3256371/files/H3K4me3_cluster_3.filtered_R1.fastq
>    https://zenodo.org/record/3256371/files/H3K4me3_cluster_3.filtered_R2.fastq
>    ```
>
> 2. Rename the files to `H3K4me3_demultiplexedR1_10000rows.fastq.gz` and `H3K4me3_demultiplexedR2_10000rows.fastq.gz`
>
> 3. Inspect the `fastq` file by using the `less` command
>

During sequencing, errors are introduced, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data.

Sequence quality control is therefore an essential first step in your analysis. We use here similar tools as described in ["Quality control" tutorial]({{site.baseurl}}/topics/sequence-analysis): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Interpretation of fastqc outputs can be found here [FastQC Output](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).



> ### Hands-on: Quality control
>
> 1. Run **FastQC** on the fastq files
>    - `fastqc H3K4me3_demultiplexedR1_10000rows.fastq.gz` and `fastqc H3K4me3_demultiplexedR2_10000rows.fastq.gz`
> 2. Inspect the generated HTML files
>
>    > ### Questions
>    >
>    > 1. How is the quality of the reads in `H3K4me3_demultiplexedR1_10000rows.fastq.gz`?
>    > 2. And in `H3K4me3_demultiplexedR2_10000rows.fastq.gz`?
>    > 3. Why do the lengths of the bases differ? Hint: Remember that the `fastq` files have been demultiplexed.
>    > 4. What should we do if the quality of the reads is not good?
>    > 5. What are the most common start sequences for the `fastq` files? How does it differ for `R1` and `R2`?
>    >
>    > > ### Solution
>    > > 1. The reads in `H3K4me3_demultiplexedR1_10000rows.fastq.gz` are of good quality
>    > >     - The mean quality score over the reads is quite high
>    > >
>    > >
>    > >
>    > >
>    > >
>    > >
>    > >
>    > >     - No N in the reads
>    > >
>    > >        ![Per base N content for read1](images/read1_per_base_n_content.png "Per base N content")
>    > >
>    > >     - No duplicated sequences
>    > >
>    > >        
>    > >
>    > >    
>    > >
>    > >        
>    > >
>    > > 2. The reads in R2 are a bit worse:
>    > >
>    > >
>    > > 3. If the quality of the reads is not good, we should:
>    > >    1. Check what is wrong and think about it: it may come from the type of sequencing or what we sequenced (high quantity of overrepresented sequences in transcriptomics data, biaised percentage of bases in HiC data)
>    > >    2. Ask the sequencing facility about it
>    > >    3. Perform some quality treatment (in a reasonable way to not loose too much information) with some trimming or removal of bad reads
>    > > 4. TODO:
>

# Step 2: Mapping of the reads

We obtain sequences corresponding to a portion of DNA linked to the histone mark of interest, H3K4me3 in this case. H3K4me3 is associated with active chromatin. It would be interesting to know if there is a difference across cell types.

{% include topics/sequence-analysis/tutorials/mapping/mapping_explanation.md
    to_identify="binding sites"
    mapper="Bowtie2"
    mapper_link="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"
    answer_3="Wang et al. (2018) did ChIP-seq on mouse. So we should use the mouse reference genome. We will use mm10 (the latest build)"
%}

## Running BWA

> ### Hands-on: Mapping
>
> 1. Run bwa on the fastq files, using an mm10 reference genome.
> TODO Add reference index files into the server.
>
> `bwa mem -t 1 $ref ${f}_R1.fastq.gz ${f}__R2.fastq.gz | samtools view -b - > $outf`
>
> 2. Inspect the mapping stats
>
>    > ### Questions
>    >
>.   > How do you view the aligned reads?
>.   > How do we see the length of each chromosome used in the mapping? Hint: see header in samtools
>    > How many reads where mapped? Uniquely or several times?
>    >

The output of BWA is a BAM file (binary of Sequence Alignment/Map).

## Inspection of a BAM file

TODO Add more things to inspect.


## Correlation between samples

We will compare genome-wide correlation of H3K4me3 and H3K4me1 for different cell clusters.

To compute the correlation between the samples we are going to to use the QC modules of [deepTools](http://deeptools.readthedocs.io/), a software package for the QC, processing and analysis of NGS data. Before computing the correlation a time consuming step is required, which is to compute the read coverage (number of unique reads mapped at a given nucleotide) over a large number of regions from each of the inputed BAM files. For this we will use the tool **multiBamSummary** {% icon tool %}. Then, we use **plotCorrelation** {% icon tool %} from deepTools to compute and visualize the sample correlation. This is a fast process that allows to quickly try different color combinations and outputs.

Since in this tutorial we are interested in assessing H3K4me3 and H3K4me1 scChIC-seq samples. To save time, we have already done that and summarized the scChIC-seq data as bigwig files, which contains the read coverage.

> ### Hands-on: Correlation between samples
>
> 1. Download the read coverage files `bigwig`.
>
> 2. Rename the files
> 3. Compare bigwigs using `multiBigwigSummary`
>    - *"Sample order matters"*: `No`
>       - {% icon param-files %} *"BAM/CRAM file"*: the 8 imported BAM files
>    - *"Choose computation mode"*: `Bins`
>       - *"Bin size in bp"*: `100000`
>
>
>    Using these parameters, the tool will take bins of 100000 bp. For each bin the overlapping reads in each sample will be computed and stored into a matrix.
>
> 4. **plotCorrelation** with the following parameters
>    - *"Correlation method"*: `Spearman`
>    - Plot `heatmap` or `scatterplot`.
>.   - Plot output in log scale (--log1p)
>
>    

> ### Questions
>
> From the correlation plot, can you infer which clusters correspond to the same cell type in H3K4me1 and H3K4me3?
>

# Step 4: Exploring `bam` and `bigwig` files on the IGV browser

The bam file contains only reads falling in specific genomic regions, in order to reduce the file size.

Cell-type specific regions to look at:

```
chr7    114972165       116898716
chr7    103325741       104425479
chr3    90523036        90870443
chr11   44114099        45269522
```



> ### Questions:
>     1. Can you infer which clusters correspond to which cell types based on the coverage around the four regions?
>

# Step 5: Detecting enriched regions (peak calling)

We could see in the scChIC-seq data some enriched regions that differ across samples. We now would like to call these regions to obtain their coordinates, using `hiddenDomains`

> ###  Hands-on: Peak calling
>
> 1. Calling peaks with `hiddenDomains`
>    > Required inputs:
>    >     -g Size of chromosomes for the mouse genome.
>    >
>    >     -q Minimum MAPQ score. What is an appropriate MAPQ score to use? A low MAPQ score may include reads that are poor quality, while high MAPQ score keeps only high quality reads. Calculate, for each `bam` file, how many reads are assigned at each MAPQ score. Hint: the MAPQ score of each read is found in the fifth column of the bam file. Hint 2: UNIX tools such as `cut` or `awk`, `sort`, `uniq` can do this calculation quickly. 
>    >
>    >     -p A threshold to remove domains called with probabilities less than `p`. Set to a value such as 0.5. You can play around with this value to see how it changes the output.
>    >     -b minimum length. Use the default 1000 bp.
>
>    > ### Questions
>    >
>    > Which type of files were generated? What do they include?
>    > How do the peaks called differ between samples? Can you infer cell types based on this analysis? 
>   


# Conclusion

We learned to explore `fastq` and `bam` files as well as do calculations on them. We visualized `bam` files and `bigwig` files to see how reads are mapped on the genome. 

