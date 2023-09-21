### Installing packages ###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# It is very important for the analysis to have the latest versions of BiocManager and DiffBind! DOUBLE CHECK!!
# However, the latest Biocmanager does not include GenomicRanges, which DiffBind requires.
# Install Genomic Ranges individually
#BiocManager::install("DiffBind")
#install.packages("rgl") # optional, for 3D illustrations

########### Load necessary libraries ###################################
library(DiffBind)
library(rtracklayer)
library(tidyverse)
library(rgl)


########### Reading input files ###################################
### reading the sample sheet ####

## generate the sample sheet with all samples, or a selection to be directly compared for differential analysis

samples <- read.csv("SampleSheet.csv")


### reading the peaks ####
ATAC_peaks <- dba(sampleSheet = "SampleSheet.csv")

ATAC_peaks ## display the metadata associated with the DBA object
# first row: total No of unique peaks, after all overlapping were merged
# "intervals": number of peaks read for each sample

## Generate sample correlation plot based on read peaks
plot(ATAC_peaks, main = "dba.read")

########### Counting reads in each interval for each sample ###################################
ATAC_count <- dba.count(ATAC_peaks, summits = 150)

ATAC_count ## display metadata
# Row one of result shows the length of consensus peakset used in all analyzed samples
# FRiP column indicates proportion of peaks per sample, which overlap a peak in the consensus set

## Generate sample correlation plots based on counted peaks in each sample
plot(ATAC_count, main = "dba.count")
dba.plotPCA(ATAC_count, attributes = DBA_ALL_ATTRIBUTES)


########### Normalisation/establishing contrast(s)/differential analysis ###################################

### applying default normalisation ###
ATAC_norm <- dba.normalize(ATAC_count)


### establishing contrasts ####
ATAC_contrast <- dba.contrast(DBA = ATAC_norm, categories = DBA_CONDITION)
ATAC_contrast


### performing differential analysis ###
bone_diff <- dba.analyze(ATAC_contrast)
plot(bone_diff, contrast=1, main = "TITLE")

## retrieving differentially open chromatin
bone_diff.DB <- dba.report(bone_diff)
bone_diff.DB


########### Plotting/data export ###################################

dba.plotMA(bone_diff, method = DBA_DESEQ2, bNormalized = TRUE, factor = "Normalized data")
dba.plotPCA(bone_diff,attributes = c(DBA_TISSUE,DBA_CONDITION),label = DBA_REPLICATE)
dba.plotPCA(bone_diff, b3D=T, attributes = c(DBA_TISSUE,DBA_CONDITION),label = DBA_REPLICATE) #3D plot incl 3rd principal component

dba.plotHeatmap(bone_diff, contrast = 1, correlations = FALSE)
dba.plotVolcano(bone_diff,contrast = 1) # volcano plot shows differentially accessible peaks


### export bed file with metadata from differential analysis, for use with HOMER ####
#access elements of the list
bone_diff_db_clean <- bone_diff.DB %>% 
  as.data.frame() %>%  # reading the metadata file as simple R dataframe
  rownames_to_column("peakID") %>% # including the Unique peak identifiers as a separate column named "peakID"
  mutate(strand=1) %>% # turn asterisks into HOMER-readable strand name
  select(seqnames,start,end,peakID,width,strand,Conc,Conc_positive,Conc_negative,Fold,p.value,FDR) #order output columns according to HOMER needs


##inspecting file
bone_diff_db_clean %>% group_by(strand) %>% summarize(count=n()) # which characters appear in strand column
bone_diff_db_clean %>% nrow  ==  bone_diff_db_clean$peakID %>% unique %>% length # is total row number of dataframe exactly equal to the number of unique peakIDs?



# HOMER fails if RefSeq chromosome names are used. Used following line to generate HOMER-readable chromosome names
tableTranslate <- read.table("refseq_to_chr.txt", sep = "", header = TRUE) #read translation table

chrom_names_all_peaks <- left_join(bone_diff_db_clean, tableTranslate, by = ("seqnames"="seqnames")) %>% # join tables according to seqnames column
  select(-"seqnames") %>% # remove seqnames column, keeping only "chromosome" column
  relocate("chromosome") # relocate "chromosome" column to become column 1

# write table containing result, including all positive and negative fold changes
chrom_names_all_peaks %>% write.table(file="UCSCbone_diff_allMetadata.bed",col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)

# separate negative from positive fold changes, write results into separate bed files
# these output files are to be used in homer analysis
bone_diff_clean_entp_neg <- chrom_names_all_peaks %>%
  filter(Fold<0) %>%
  write.table(file = "entpd_neg_enriched.bed",col.names = FALSE, sep = "\t",quote = FALSE, row.names = FALSE)

bone_diff_clean_entp_pos <- chrom_names_all_peaks %>%
  filter(Fold>0) %>%
  write.table(file = "entpd_POS_enriched.bed",col.names = FALSE, sep = "\t",quote = FALSE, row.names = FALSE)

