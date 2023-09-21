#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

##Once the core Bioconductor packages are installed, we can install the topGO package by

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler",force = T)
#BiocManager::install("HDO.db",force = T)
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Dr.eg.db")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Dr.eg.db)
library(ggplot2)

##loading significantly regulated genes

sig <- read.csv("SIG.csv", row.names = 1)

##selecting only genes with >0.5 log2 fold change (enriched in entpd5a+ population)

sig <- sig[sig$log2FoldChange > 0.5,]

## getting row names (Ensembl IDs of corresponding genes)

genes_for_GO <- rownames(sig)

##run GO

GO_res <- enrichGO(gene = genes_for_GO, 
                        OrgDb = "org.Dr.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")

## convert to dataframe & plot

GO <- as.data.frame(GO_res)
write.csv(GO,file = "GO.csv",row.names = FALSE)

ggplot(GO[1:20,], aes(x = reorder(Description, Count), 
                           y = Count, 
                           fill = p.adjust))+
  geom_bar(stat = "identity", color="black", width=0.4)+
  theme_bw()+
  ylab("")+ 
  xlab("")+
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.text.x = element_text(size=10)) + # setting font size
  coord_flip()


            