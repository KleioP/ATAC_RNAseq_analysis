##load libraries
library(tidyverse)
library(plotly)
library(VennDiagram)
library(RColorBrewer)

#########################  Pre-processing Zv11 gene dataset for integrative analysis ####################### 

# read Zv11 assembly: all transcripts with their locations
# mart_export_Zv11.csv is generated using the Ensembl Biomart tool
# optional: remove different transcript entries, when gene location is the same
all_genes <- read.csv(file = "mart_export_Zv11.csv", header = T, sep = ",")
all_genes_merge <- all_genes %>%
  select(Transcript.ID, Gene.name, Chromosome,
         Transcript.start, Transcript.end, Strand) #%>%
#distinct() #Keep only unique/distinct rows # optional

#isolating all genes in the +1 strand in a file
all_genes_plus <- all_genes_merge %>% filter(Strand == 1)
#isolating all genes in the -1 strand in a file
all_genes_minus <- all_genes_merge %>% filter(Strand == -1)

# for the + strand genes:
# set promoter.start column: required distance before transcript start
# set promoter.end column: end of transcript. This is to cover potential intron regulatory regions

distance_from_promoter <- 5000 # (bp) set this parameter according to requirement, here 5kb

all_plus_prom <- all_genes_plus %>%
  mutate(promoter.start = Transcript.start - distance_from_promoter) %>%
  mutate(promoter.end = Transcript.end) %>%
  select(Transcript.ID, Gene.name, Chromosome,
         peak.region.start = promoter.start,
         peak.region.end = promoter.end,
         TSS=Transcript.start, Strand)
  #using select() to reorder and rename columns at the same time

# similar processing for - strand genes
# taking account of the fact that Transcript.Start here indicates the end of the 3' UTR
# and Transcript.end indicates the TSS

all_minus_prom <- all_genes_minus %>%
  mutate(promoter.start = Transcript.end + distance_from_promoter) %>%
  mutate(promoter.end = Transcript.start) %>%
  select(Transcript.ID, Gene.name, Chromosome,
         peak.region.start = promoter.end,
         peak.region.end = promoter.start,
         TSS = Transcript.end, Strand)

#merge vertically plus and minus promoters
mergepromoters <- bind_rows(all_plus_prom,all_minus_prom)
write.csv(mergepromoters, file = "5kb_upstream_all_genes_promoter_areas.csv",
          row.names = F)


#########################  Pre-processing RNAseq DEG datasets for integrative analysis ####################### 

# reading in datasets
DEGs <- read.csv(file = "{rnaseq_results.csv}", header = T)

# filtering positive vs negative fold changes
# incorporating gene names to the tables with left_join()
# keeping only the unique entries for each gene with distinct()

RNAseq_pos <- DEGs %>% 
  filter(log2FoldChange > 0.5) %>%
  rename(Gene.ID=X) %>%
  left_join(y = all_genes, by = "Gene.ID") %>%
  select(Gene.ID, baseMean, log2FoldChange, padj, Gene.name) %>%
  distinct()
write.csv(pos, "pos.csv", row.names = F)

RNAseq_neg <- DEGs %>% 
  filter(log2FoldChange < -0.5) %>%
  rename(Gene.ID=X) %>%
  left_join(y = all_genes, by = "Gene.ID") %>%
  select(Gene.ID, baseMean, log2FoldChange, padj, Gene.name) %>%
  distinct()
write.csv(neg, "neg.csv", row.names = F)

#########################  Intersecting RNAseq datasets from different cell types ####################### 

# run first the script Functions.R
## Run "Process_Venn_Data" custom function to extract gene names column as single vector
# for details on function "Process_Venn_Data" see script "Functions.R"

RNAseq_pos_gene_names <- Process_Venn_Data(RNAseq_pos)
RNAseq_neg_gene_names <- Process_Venn_Data(RNAseq_neg)
all <- Process_Venn_Data(all_genes)

# plotting Venn Diagram comparing positive and negative regulated genes
 # in comparison to all genes in assembly
# Similarly, outputs from different samples can be compared

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x=list(RNAseq_neg_gene_names, RNAseq_pos_gene_names, all),
  category.names = c("Neg", "Pos", "All"),
  filename = "{file_name.tiff}",
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
