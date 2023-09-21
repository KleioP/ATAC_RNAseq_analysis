# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(org.Dr.eg.db)
library(AnnotationDbi)

## load differentially expressed gene (DEG) list

sig <- read.csv("All_posit_vs_negative_pAdj005.csv", row.names = 1)
sig_df <- as.data.frame(sig)

## add a "symbol" column on the dataframe containing the gene symbol for each ENS ID

sig_df$symbol <- mapIds(org.Dr.eg.db, keys = rownames(sig_df), keytype = "ENSEMBL", column = "SYMBOL")

## filter rows that will appear in heatmap
# ensuring padjusted < 0.05 has been filtered
# baseMean of all samples for that gene > 50 (to avoid noise)
# filter for higher logfold changes by setting the absolute value > 2
sig_filt <- sig_df[(sig_df$padj < 0.05) & 
                     (sig_df$baseMean > 50) & 
                     (abs(sig_df$log2FoldChange) > 2),] 

# order by logfold2 change
sig_filt <- sig_filt[order(sig_filt$log2FoldChange, decreasing = T),]

# use dds object from DESEQ2
# get normalised counts from dds object
# this is also done in the DESEQ2 script,  see vignette

rld <- rlog(dds, blind = F)

# make matrix from rld object
# only contains rows present in sig_filt
# samples_sheet: object created in DESEQ2 analysis, contains sample names
mat <- assay(rld)[rownames(sig_filt), rownames(sample_sheet)] 

colnames(mat) <- rownames(sample_sheet)
base_mean <- rowMeans(mat)

# center and scale each column using scale() function
# aka. get Z-scores
# the transpose matrix using t() function
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- colnames(mat)

# set number of genes to keep from sig_filt dataframe for plotting 
num_keep <- 25
# like this we create a sequence which is a combination of sequences, to take 
  # the top n=num_keep rows from mat.scaled and the last values of the mat.scaled
  # starting n=num_keep before the end of the total rows, until the end of the rows
rows_keep <- c( seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

# take up log2 values for the genes we are keeping (only from chosen rows)
log2_val <- as.matrix(sig_filt[rows_keep,]$log2FoldChange)
colnames(log2_val) <- "log2FC"

# similarly, get corresponding means for each gene we are keeping
mean <- as.matrix(sig_filt[rows_keep,]$baseMean)
colnames(mean) <- "AverExp"

# make colour map for logfold expression
# maps values between 3 colours for min and max log2 values
col_logFC <- colorRamp2(c(min(log2_val), 0, max(log2_val)), c("blue", "white", "red"))

# make a map for average expression between white and red
  # maps between 4 quantiles of mean values:
  # quantile(mean) returns 4 quantiles for the values of mean
  # choosing quntile(mean)[1] returns the first quantile etc
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

## make heatmap

# heatmap annotations
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2),
                                               height = unit(2, "cm")))

# three components of heatmap h1, h2, h3

# h1: Zscores

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F,
              column_labels = colnames(mat.scaled), name = "Z-score")

# h2: log Fold Change

h2 <- Heatmap(log2_val, 
              #row_labels = rownames(sig_filt[rows_keep,]), # for ENSEMBL symbols to show
              row_labels = sig_filt$symbol[rows_keep], # altern, for gene symbols to show
              cluster_rows = F, name = "logFC",
              #top_annotation = ha, # optional if we want top box to show
              col = col_logFC, # assigning colour map
              cell_fun = function(j, i, x, y, w, h, col) {  # adding labels, optional step
                grid.text(round(log2_val[i, j], 2), x, y)
              })

# h3: Average expression

h3 <- Heatmap(mean,
              #row_labels = rownames(sig_filt[rows_keep,]), # for ENSEMBL symbols to show
              row_labels = sig_filt$symbol[rows_keep], # altern, for gene symbols to show,
              cluster_rows = F,
              name = "AverExp",
              col = col_AveExpr, # assigning colour map
              #cell_fun = function(j, i, x, y, w, h, col) {  # adding labels, optional step
              #  grid.text(round(mean[i, j], 2), x, y)}
              ) # adding labels, optional step
              
# combine heatmaps

h <- h1 + h2 + h3


# saving object

tiff("./heatmap_all_880genes.tiff", res = 600, width = 5000, height = 10000)
print(h)
dev.off()
