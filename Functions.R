## the following function extracts a column y (default: "Gene.name") from a datafrane x

Process_Venn_Data <- function(x, y = "Gene.name"){
  # function with arguments a given dataframe (x), and a column of that dataframe (y)
  
 #The extracted vector is made up only of unique values present in given column 
  vector <- unique(x[[y]])
  
  # Converting any missing values in the Vector to NAs
  vector[vector == ""] <- NA
  
  # filtering for NA and keeping only non NAs
  filt <- !is.na(vector)
  clean_vector <- vector[filt]
  
  return(clean_vector)
}


## the following function is used to extract genes that contain ATAC peaks within 
   # given peak.region.start/end

Genes_with_Peaks <- function(x, y){
  #Arguments: 
    #x: dataframe with genes of interest and their associated peak.start
      # and peak.end
    #y: processed enrichment ATAC file for a cell population (renamed, 
      # mutated, selected columns)
  
  #full_joining the promoter.region.start/end dataset with the peak enrichment 
    # file by Chromosome
  all_prom_all_peaks <-full_join(x, y , "Chromosome")
  
  # filtering for only the rows where ATAC "peak.start/end" fall within the defined 
  # "peak.region.start/end" of each transcript:
  # a. peak.region.start < peak.start < peak.region.end
  # b. peak.region.start < peak.end < peak.region.end
  # Both (a) & (b) must be true
  genes_with_peaks <- all_prom_all_peaks %>%
    rowwise() %>%
    filter(between(peak.start, peak.region.start, peak.region.end) & 
             between(peak.end, peak.region.start, peak.region.end))
  
  return(genes_with_peaks)
  
}
