setwd("/Users/mariapau/desktop/master/bioinfo/1st_TERM/PGB/bHLH_project/ChIP_seq/IP_meme")

# CREATING TABLES FOR:
## Motif enrichment

library(readr)
library(dplyr)
library(stringr)

?write.table()

summary <- read.delim("summary.tsv", sep="\t", header=TRUE) #Read the summary.tsv
summary_simple <- summary [, c("MOTIF_SOURCE", "MOTIF_ID", "ALT_ID", "CONSENSUS", "SITES", "E.VALUE", "E.VALUE_SOURCE")]
write.table(summary_simple, "summary_simple.tsv",sep="\t",row.names=FALSE,quote=FALSE)
View(summary_simple)



## Motif centrality

summary2 <- read.delim("centrimo_out/centrimo.tsv", sep="\t", header=TRUE, check.names = FALSE) #Read the summary.tsv
summary_simple2 <- summary2 [, c("motif_id", "motif_alt_id", "consensus", "E-value", "bin_location", "sites_in_bin", "total_sites", "p_success")]
write.table(summary_simple2, "summary_simple2.tsv",sep="\t",row.names=FALSE,quote=FALSE)

View(summary_simple2)

# CHIPSEEKER: All the steps I have used to do it are in https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

#### installing packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install("clusterProfiler")

BiocManager::available() #To see which packages are there

BiocManager::install("ChIPseeker")

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

BiocManager::install("org.Mm.eg.db")

BiocManager::install("ReactomePA")

#### loading packages

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene #contains all the gene coordenates (annotation database). Without this, ChIPseeker would not know where genes or promoters are. 
library(clusterProfiler)

## ---------- CHIPSEEKR ANALYSIS -----------
setwd("/Users/mariapau/desktop/master/bioinfo/1st_TERM/PGB/bHLH_project/ChIP_seq")

### ChIP profiling --------
peak <- readPeakFile("IP_summits.bed") #readPeakFile loads the peak and store in GRanges object
peak

### ChIP peaks coverage plot --------
covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))


### Profile of ChIP peaks binding to TSS regions ---------
#### Prepare the TSS regions and then align the peaks that are mapping to these regions, and generate the tagMatrix. 
#### to speed up the compilation of this vignettes, we use a precalculated tagMatrix

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tm <- getTagMatrix(peak, windows = prom)#tagMatrix
tagHeatmap(tm)  # heatmap only

?peakHeatmap()
peakHeatmap(
  peak = peak,
  TxDb = txdb,
  upstream = 3000,
  downstream = 3000,
  by = "gene",
  type = "start_site",
  nbin = 20  # 200–800 is typical; reduce if plots are too smooth/slow
)

sum(overlapsAny(peak, genes(TxDb.Mmusculus.UCSC.mm10.knownGene)))#Test how many peaks overlap mm10 gene regions  


### Average Profile of ChIP peaks binding to TSS region ---------
plotAvgProf(tm, xlim = c(-3000, 3000), conf = 0.95, resample = 1000) 

plotAvgProf2(
  peak = peak, 
  TxDb=txdb, 
  upstream=3000, 
  downstream=3000,
  xlab="Genomic Region (5'->3')", 
  ylab = "Read Count Frequency",
  conf = 0.95
  )

## Peak Annotation -------------
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)  # mm10
library (org.Mm.eg.db)                         # mouse gene annotations
library(GenomeInfoDb)
peakAnno <- annotatePeak(
  peak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mm.eg.db"
)
plotAnnoPie(peakAnno)
vennpie(peakAnno)
#install.packages("ggupset")
#install.packages("ggimage")
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


## Functional enrichemnt analysis ------
library(ReactomePA)
library(org.Mm.eg.db)
library(tidyr)
### Extract Entrez IDs from annotatePeak result
gene_ids <- as.data.frame(peakAnno)$geneId
gene_ids <- unique(na.omit(gene_ids))

### Reactome enrichment (mouse)
path_mouse <- enrichPathway(gene = gene_ids,
                            organism = "mouse",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                            readable = TRUE)   # converts Entrez -> SYMBOL using org.Mm.eg.db

head(path_mouse, 3)
#### Plots
barplot(path_mouse, showCategory = 20)
dotplot(path_mouse, showCategory = 20)


## Extract gene IDs involved in the enriched pathway :)
#slotNames(path_mouse)
enriched_genes <- path_mouse@result #Extracting the "result" slot inside that object
head(enriched_genes)

### Separate the gene IDs for each enriched pathway
enriched_genes_df <- enriched_genes %>%
  select(ID, Description, geneID) %>%
  mutate(geneID = strsplit(geneID, "/")) %>%  # Split gene IDs by "/"
  unnest(geneID)  # Flatten the list into a dataframe

#### Display the resulting dataframe
head(enriched_genes_df)

#### Save the gene IDs to a CSV file
write.csv(enriched_genes_df, "enriched_genes_mouse.csv", row.names=FALSE)

# Get Ensembl IDs ----------
## By extractin the Esembl IDs of the found genes, they could be compared with the RNA-seq analysis. 
library(org.Mm.eg.db)
library(dplyr)

## Creation of a list contianing Ensembl gene identifiers
#head(enriched_genes_df$geneID) # to check what is the IDs really look like

entrez_ids <- enriched_genes_df$geneID
ensembl_ids <- mapIds(
  org.Mm.eg.db,
  keys = entrez_ids,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first" #If multiple mappings, take the first one!
)

## Add the Ensembl IDs as a new column to the existing dataframe
enriched_genes_df$Ensembl_ID <- ensembl_ids

## Display the updated dataframe
head(enriched_genes_df)

write.csv(enriched_genes_df, "enriched_genes_with_ensembl.csv", row.names=FALSE)




