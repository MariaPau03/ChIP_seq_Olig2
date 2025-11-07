# ChIP-seq analysis
The ChIP-seq workflow integrates experimental immunoprecipitation with computational analysis to identify DNA regions bound by specific proteins. In this study, DNA fragments associated with the Olig2 mouse transcription factor were isolated through chromatin immunoprecipitation and sequenced using Ilumina technology. The resulting reads were subjected to a standard bioinformatics pipeline, which included quality control with FastQC, alignment to the Mus musuclus reference genome (mm10) using BWA and conversion of alignment files to stored and indexed BAM format with Samtools. Enriched binding regions (peaks) were then identified with MACS2 using the input control sample for background correction, followed by peak annotation and visualization to interpret the genomic distribution of bHLH binding sites. 
To start with, publicly available ChIP-seq datasets corresponding to bHLH transcription factor binding in Mus musculus were obtained from the Sequence Read Archive (SRA). Two datasets were analyzed, SRR1583894 and SRR1583890. By downloading the SRA toolkit we obtained two fastq files with only 1M reads (in order to not spend too much time with the whole total reads). 

```
$ fastq-dump -X 1000000 --split-files -v --gzip --outdir ./ SRR1583890
$ fastq-dump -X 1000000 --split-files -v --gzip --outdir ./ SRR1583894
```

By typing this part of code, it is possible to see how many reads each sample has, but as explained previously, each sample has around 1M reads, as detailed in the previous code. 

```
$ gzcat SRR1583894_1.fastq.gz | grep -c "^@" 
$ gzcat SRR1583890_1.fastq.gz | grep -c "^@" 
```
Continuing with the quality control, raw sequencing data were evaluated using FastQC to assess per-base quality, GC content, sequence duplication and adapter contamination. The following code is necessary to have FastQC executable and add it in the variable $PATH.

```
$ wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip 
$ unzip fastqc_v0.12.1.zip
$ chmod 770 -R FastQC
$ export PATH=/Applications/Softwares/FastQC:$PATH
$ fastqc SRR1583894_1.fastq.gz
$ fastqc SRR1583890_1.fastq.gz
```
Among the output files, one of the most important things was the html file, in which BoxWhisker plots were represented (as shown in Figure x). 
Secondly, all reads were aligned to the mouse reference genome (mm10), downloaded from the UCSC Genome Browser. Furthermore, this is one of the core steps in the ChIP-seq analysis pipeline. Before alignment, the reference was indexed with BWA. 

`$ bwa index mm10.fa`

Sequencing reads were mapped to the reference genome using BWA-MEM. For single-end data, the following command was used :

```
$ bwa mem mm10.fa  SRR1583890_1.fastq.gz > SRR1583890_1.sam
$ bwa mem mm10.fa  SRR1583894_1.fastq.gz > SRR1583894_1.sam
```
Take into account that his process might take some minutes to run. What it basically does is using the mem algorithm which will find the maximal exact match for short-read alignment and single-end data, in other words, it finds where each sequencing read matches best in the genome. The output will be redirected into a .sam file which will contain its name, where it aligns in the genome, alignment score, mapping quality and CIGAR string (insertions, deletions and matches). Without this step, we wouldn’t have been able to find where in the genome each read was originated from, and all this information is the foundation for the later peak calling, visualization and motif discovery among other things.

```
$ samtools sort -O bam SRR1583890_1.sam > SRR1583890_1.bam
$ samtools sort -O bam SRR1583894_1.sam > SRR1583894_1.bam
```
The following code is necessary to create .bam.bai files in order to visualize the bam files in the IGV viewer.

```
$ samtools index SRR1583890_1.bam 
$ samtools index SRR1583894_1.bam  
```
Peaks representing genomic regions enriched in bHLH binding were identified using MACS2. The ChIP sample was compared against its input control to account for background signal.

```
$ pip install MACS2 
$ macs2 callpeak -g mm -f BAM -t SRR1583890_1.bam -c SRR1583894_1.bam --bw 200 --outdir . -n IP`
```

Peak calling produced a list of significant binding regions in BED format, which were later used for visualization and annotation. What is more, it has been done an intersection between the human Olig2 to find the overlapping peaks between the two experiments, in other words, to get more robust peaks, more reliable results. For that, it was used Bedtools to find overlapping peaks between the two experiments.  

`$ intersectBed -a IP.peaks.narrowPeak -b IP_2.narrowPeak >RobustPeaks.bed`

Then, to know which motifs are centrally enriched with the peaks, it is needed a previous step to resize peaks to equal length form the midpoint and once they have the same size and have obtained the sequences from the reference genome using again Bedtools, install meme-chip as shown.

```
$ awk '{print $1,$2+int(($3-$2)/2)-250,$2+int(($3-$2)/2)+250}' RobustPeaks.bed | tr " " "\t" | sort -k1,1 -k2,2n > RobustPeaks.resized.bed
$ fastaFromBed -fi mm10.fa -bed RobustPeaks.resized.bed > RobustPeaks.resized.bed.fa
$ conda install -c bioconda meme
```

To run meme-chip, is necessary to have the jaspar.meme which is a file containing known motifs from the JASPAR database and the Robust.Peaks.resized.bed.fa (input FASTA file with the sequences from the peaks around each summit midpoint. 

`$ meme-chip -oc IP_meme --db jaspar.meme RobustPeaks.resized.bed.fa`

To finish with, the identified peaks were annotated to the nearest genes using ChIPseeker (R package) and visualized in IGV. Coverage profiles and peak distributions were also examined using deepTools to asses genome-wide enrichment patterns. 

Last but not least, the ChIPSeeker analysis with R/Bioconductor is also attached (ChIP_seq.r) and all the results are inside the ChIP_seeker folder. 

PD: There's also two csv files in which shows enrichment output into a tidy, analysis-ready format and added Ensembl IDs to make it robust for any downstream integration or visualization.
