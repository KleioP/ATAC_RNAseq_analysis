#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")

########### Load necessary libraries ###################################
library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)


########### Preparation for diff. expression analysis ###################################

### Generating a matrix with genes and counts for all samples ### 
#TODO THese should better be in subfolder
all_files <- list.files(pattern = "*.txt") #filtering feature count output files present in working directory

list <- list() # generating empty list

for(i in all_files){
  list[[i]] <- read.csv(file = i, skip = 1, header = TRUE, sep = "\t")
  # writing content of all counts files into the list
}

# deleting unwanted columns from all list elements, keeping only "Geneid" and "counts"
list_select_col <- lapply(list, function(x) x[!names(x) %in% c("Chr","Start","End","Strand","Length")])

# use reduce function to generate dataframe by combining the Geneid columns of all elements of the list
# full_join function used within the reduce function to ensure keeping of all gene names present in each element of the list
# piping into column_to_rownames to turn column "Geneid" into rownames to allow matrix creation
x <- purrr::reduce(list_select_col, full_join, by=c("Geneid"="Geneid")) %>% 
  column_to_rownames(var = "Geneid")

# converting dataframe to matrix to use in DESeq2
RNAseq_counts <- as.matrix(x)
# turning all colnames into the all_files names
colnames(RNAseq_counts) <- all_files

#TODO ideally not hard coded. looks more profressional
#reordering columns to allow for subsetting by comparison in later step
col.order <- c("cartilRep1.txt", "cartilRep2.txt", "cartilRep3.txt",
               "osteobRep1.txt", "osteobRep2.txt", "osteobRep3.txt", 
               "intersegRep1.txt", "intersegRep2.txt","intersegRep3.txt",
               "segmRep1.txt", "segmRep2.txt", "segmRep3.txt")
RNAseq_counts_order <- RNAseq_counts[,col.order]



########### Differential Expression Analysis ###################################

## Comparison between osteoblast and cartilage replicates

# if individual tissues are to be analysed
# subsetting RNAseq_counts_order
#head_posit_vs_negative <- RNAseq_counts_order[,1:6]
#trunk_posit_vs_negative <- RNAseq_counts_order[,7:12]

#reading sample sheet
sample_sheet <- read.csv("coldata_sample_sheet.csv", row.names = 1)

# create deseq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = RNAseq_counts_order,
                              colData = sample_sheet,
                              design = ~ condition)

# prefiltering of low reads
# removal of all rows with sum of gene counts accross row < 10
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]


#setting levels for comparison
dds$condition <- factor(dds$condition, levels = c("negative","positive"))

#differential expression analysis
dds <- DESeq(dds)
res <-results(dds, alpha = 0.05) #adjusting FDR filtering cutoff for pvalue to 0.05 from 0.1
write.csv(res,file = "All_posit_vs_negative.csv")
plotMA(res,ylim=c(-2,2), main="All_posit_vs_negative")

# filtering out padj values >0.05 for later use
filt <- read.csv("All_posit_vs_negative.csv")
keep_filt <- filt[,7]<=0.05
res_p005 <- filt[keep_filt,]
filt_noNA <- na.omit(res_p005) #omitting rows with padj = NA (due to independent filtering of results function)
write.csv(filt_noNA,file = "All_posit_vs_negative_pAdj005.csv", row.names = FALSE)


####optional step####
# retrieving normalised counts for later use 
# (step included in DESeq function later, not needed for DE analysis)
dds_norm <- estimateSizeFactors(dds) #adds size factors info for each analysed sample
sizeFactors(dds_norm) #view size factors
norm_counts <- counts(dds_norm, normalized=TRUE) #retireving normalized counts
write.table(norm_counts, file="norm_counts.txt", sep="\t", quote=F, col.names=NA)



########### Plotting individual genes ###################################


# plot individual genes
# see below for automated plotting
plotCounts(dds, gene="ENSDARG00000019516",
           intgroup="condition", main = "osterix/sp7",
           xlab = "entpd5a expression")


### exploratory automated plotting of individual gene counts in both tissues
# the file genes_for_plotting.csv contains Ensembl IDs of genes of interest
genes_to_plot <- read.csv("genes_for_plotting.csv")
for(i in 1:nrow(genes_to_plot)){
  plotCounts(dds, gene=genes_to_plot[i,1],
             intgroup="condition", main = genes_to_plot[i,2],
             xlab = "entpd5a expression")
}

# plotting genes of interest with ggplot for better images
osx <-plotCounts(dds, gene="ENSDARG00000019516", 
                 intgroup="condition", returnData = TRUE)

ggplot(osx, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400,1500)) +
  xlab("cells sorted for entpd5a") +
  ylab("Gene count")+
  theme_classic()


########### Adjustment for MA plots ###################################


### Log fold change shrinkage for visualization and ranking
# the large fold changes from genes with lots of statistical information are 
# not shrunk, while the imprecise fold changes are shrunk
## following DESeq2 vignette
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_positive_vs_negative", type="apeglm")
resLFC
plotMA(res,ylim=c(-2,2),main="no adjustment")
plotMA(resLFC,ylim=c(-2,2),main="LFC")



########### Plotting Heatmaps/PCA to QC samples ###################################
### here heatmaps are for assessment of correlation between samples

### Value transformation prior to heatmap plotting
# this is done to account for inherently significant variations around the mean
#   observed when counts are low
# here 2 methods are tested: vsd = variance stabilizing transformation and 
#   rlog = regularised log and compared to a basic transformation (made using
#   normTransform)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3) # visualisation


# the following gives log2(n + 1), where n=gene count
ntd <- normTransform(dds) # creates a DESeqTransform object

# in each case, visualising standard dev across counts
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


# generating heatmaps using all 3 methods of value transformation
select <- order(rowMeans(counts(dds,normalized=TRUE)), # select top genes
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"condition"])
rownames(df) <- colnames(assay(ntd)) #give rownames to df (for pheatmap to work)

pheatmap::pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_col=df)
pheatmap::pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_col=df)
pheatmap::pheatmap(assay(rld)[select,], cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_col=df)

## heatmap of sample correlation
#getting sample distances
sampleDists <- dist(t(assay(rld))) # apply the dist function to the transpose of the transformed count matrix 

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(assay(rld))
#rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(c("blue", "white", "red"))(200)
heatmap <- pheatmap::pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="")

tiff(filename = "./allSample_corr_heatmap.tiff", res = 2000, width = 15000, height = 10000)
print(heatmap)
dev.off()


###PCA

PCA_classic <-plotPCA(rld, intgroup="condition") +
  #geom_text(aes(label=name),vjust=2,check_overlap = F,size = 2) + #run this for labels
  theme_classic() + 
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 15))
  
  
tiff(filename = "./PCA_classic.tiff", res = 2000, width = 10000, height = 10000)
print(PCA_classic)
dev.off()
