##load libraries
library(tidyverse)
library(plotly)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)

# This script first identifies ATACseq peaks for sample of interest 
 # that are found in
  # (a)  promoters of all genes in Zv11
  # (b)  promoters of the DEGs for the sample of interest

# First make sure to run the scripts "Integr_ATAC_RNAseq_DataPreparation.R"
 # and "Functions.R"

#########################  Comparing ATACseq to full Zv11 ####################### 

# reading in files DiffBind output files (.bed)
# renaming columns for coordinates
# removing "chr" from all rows in the chromosome column,
 # to make them compatible with all_genes
# keeping needed columns
ATAC_peaks <- read.csv(file = "{Peak_file.bed}", header = F, sep="\t")
ATAC_peaks <- ATAC_peaks %>%
  select(Chromosome_chr = V1, peak.start = V2, peak.end = V3)%>%
  mutate(Chromosome = str_sub(Chromosome_chr, 4, -1)) %>%
  select(Chromosome, peak.start, peak.end)


# reading Zv11 peak.region.start/peak.region.end file(s) produced in
 # "Integr_ATAC_RNAseq_DataPreparation.R"
# here, looking for peaks within 5kb, 10kb or 20kb of Transcription Start Site

all_prom_5kb <- read.csv("5kb_upstream_all_genes_promoter_areas.csv", header = T)
all_prom_10kb <- read.csv("10kb_upstream_all_genes_promoter_areas.csv", header = T)
all_prom_20kb <- read.csv("20kb_upstream_all_genes_promoter_areas.csv", header = T)


## Extracting the genes that contain ATAC peaks within
  # 5kb, 10kb and 20kb of the TSS, or in their introns
  # for details on function "Genes_with_Peaks" see script "Functions.R"

allGenes_peaks_5kb <- Genes_with_Peaks(all_prom_5kb, ATAC_peaks)
allGenes_peaks_10kb <- Genes_with_Peaks(all_prom_10kb, ATAC_peaks)
allGenes_peaks_20kb <- Genes_with_Peaks(all_prom_20kb, ATAC_peaks)

write.csv(allGenes_peaks_5kb, "Intersect_Zv11_ATAC_5kb.csv", row.names = F)
write.csv(allGenes_peaks_10kb, "Intersect_Zv11_ATAC_10kb.csv", row.names = F)
write.csv(allGenes_peaks_20kb, "Intersect_Zv11_ATAC_20kb.csv", row.names = F)

## plotting results

# extract columns with transcript IDs as single vectors
# for details on function "Process_Venn_Data" see script "Functions.R"

Intersect_ATAC_5kb <- Process_Venn_Data(allGenes_peaks_5kb, "Gene.name")
Intersect_ATAC_10kb <- Process_Venn_Data(allGenes_peaks_10kb, "Gene.name")
Intersect_ATAC_20kb <- Process_Venn_Data(allGenes_peaks_20kb, "Gene.name")

## producing bar chart indicating proportion of Zv11 genes associated with ATAC peaks
# how many have peaks found within 5kb of TSS+introns
# vs how many have peaks within 5-10kb of TSS+introns
# vs how many have peaks within 10-20kb of TSS+introns

## Reading Proportion_ATAC_coverage_on_allGenes.csv file containing the following numbers:
  # % genes with peaks within 5kb+introns out of genes with peaks within 20kb of TSS+introns
  # % genes with peaks within 5kb to 10kb+introns out of genes with peaks within 20kb of TSS+introns
  # % genes with peaks within 10kb to 20kb+introns out of genes with peaks within 20kb of TSS+introns
## calculations not shown 

Allgenes_coverage <- read.csv("Proportion_ATAC_coverage_on_allGenes.csv")

plot_allG_coverage <- ggplot(Allgenes_coverage,
                            aes(x=Cell_Type, y = Proportions,
                                fill = Group, label = Proportions))+
  # y= Number_of_Genes if actual numbers wanted
  geom_bar(stat = "identity")+
  theme_bw()+
  ylab("% of genes with ATAC peaks")+
  xlab("")+
  geom_text(size=2, position = position_stack(vjust = 0.5))

plot_allG_coverage

tiff(filename = "Proportion_allGene_coverage.tiff",width = 3000, height = 1500,res = 600)
print(plot_allG_coverage)
dev.off()

# generate Venn Diagram showing overlap between DEGs and 
 # Zv11 genes associated with ATAC peaks within 20kb of TSS+introns 

myCol <- c("red", "blue")
venn.diagram(
  x=list(Intersect_ATAC_20kb, RNAseq_pos),
  category.names = c("ATAC", "RNAseq"),
  filename = "RNAseq_vs_ATAC_segments.tiff",
  output = TRUE,
  
  # Output features
  imagetype="tiff",
  height = 1200, 
  width = 1200, 
  resolution = 400,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold"
)

#########################  Comparing ATACseq to DEGs for a sample ####################### 


# using x= RNAseq datasets (preparation Script) and
 # y = allGenes_peaks_xkb files generated above
 # going to inner_join(x,y, by = Gene.name),
 # keeping as a result only DEG entries that have ATAC peaks
 # here done for peaks within 5kb of promoter

DEGs_with_peaks <- inner_join(RNAseq_pos, allGenes_peaks_5kb, "Gene.name")

## Extract relevant peaks (located within 5kb of DEGs' TSS regions)
# adjust to HOMER readable format

Peaks_in_DEGs <- DEGs_with_peaks %>%
  select(Chromosome, peak.start, peak.end, Strand) %>%
  unique()
Peaks_in_DEGs$Chromosome <- sub("^", "chr", Peaks_in_DEGs$Chromosome)


# write bed files for use with HOMER
write.table(Peaks_in_DEGs, "Peaks_in_DEGs_for_HOMER.bed", quote = F, sep="\t", col.names = F)

# Isolate gene names in each case
# Unique names
# Remove missing values

Gene_names <- Process_Venn_Data(DEGs_with_peaks)

# saving as .csv files the list of DEGs that have at least one ATAC peak within
  # 5kb of the TSS and their introns

write.csv(Gene_names, file = "DEGs_with_peaks.csv", row.names = F)

# saving the environment for later
save.image(file='ATAC_integration_environment.RData')
load('ATAC_integration_environment.RData')

## similar preparation of RNAseq DEG datasets
## Numbers produced in script "Integr_ATAC_RNAseq_DataPreparation.R"

# plotting
## ggplot_gene coverage

## Reading Proportion_ATAC_coverage_on_DEGs.csv file containing the following numbers:
  # % genes with peaks within 5kb+introns out of DEGs with peaks within 20kb of TSS+introns
  # % genes with peaks within 5kb to 10kb+introns out of DEGs with peaks within 20kb of TSS+introns
  # % genes with peaks within 10kb to 20kb+introns out of DEGs with peaks within 20kb of TSS+introns
## calculations not shown 

DEG_coverage <- read.csv("Proportion_ATAC_coverage_on_DEGs.csv")

plot_DEG_coverage <- ggplot(DEG_coverage,
                            aes(x=Cell_Type, y = Proportions,
                                fill = Group, label = Proportions))+
                                # y= Number_of_Genes if actual numbers wanted
  geom_bar(stat = "identity")+
  theme_bw()+
  ylab("% of DEGs with ATAC peaks")+
  xlab("")+
  geom_text(size=2, position = position_stack(vjust = 0.5))

plot_DEG_coverage

tiff(filename = "Proportion_DEG_coverage_noText.tiff", width = 3000, height = 1500, res = 600)
print(plot_DEG_coverage)
dev.off()
