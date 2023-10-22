######################### ONE OFF PREP MITOCHONDRIA SEQ FILE #########################


#convert mitoch reads (fasta format) to .bam by aligning to danio rerio genome
#-f specifies fasta input
#-U for unpaired data to be aligned
bowtie2 -f -x genome -U Danio_mitoch.fasta -S Danio_mitoch.sam
samtools view -Sb Danio_mitoch.sam > Danio_mitoch.bam
samtools sort Danio_mitoch.bam > Danio_mitoch_sorted.bam
samtools index Danio_mitoch_sorted.bam > Danio_mitoch_sorted.bam.bai


######################### Loading Dependencies & Packages #########################
module add foss
module add palma/2019a
module add icc/2019.1.144-GCC-8.2.0-2.31.1
module add GCC/8.2.0-2.31.1
module add OpenMPI/3.1.3
module add impi/2018.4.274
module add ifort/2019.1.144-GCC-8.2.0-2.31.1

module add bcl2fastq2/2.20.0
module add fastp
module add Bowtie2/2.3.5.1
module add SAMtools/1.9
module add deepTools/3.3.1-Python-3.7.2
module add BEDTools/2.28.0
module add Sambamba/0.7.1
module add MACS2/2.2.6-Python-3.7.2

## to run ataqv
# module add palma/2021a
# module add GCC/10.3.0
module add ataqv/1.3.0

######################### Demultiplexing #########################

# generate .csv file (i.e. Sample Sheet) using the Illumina Experiment Manager program
# save the Sample Sheet in the first directory of the sequencing folder
# use following command to concatenate .bcl files from different sequencing lanes, and sorts by given adapter sequence

bcl2fastq -R [SEQUENCING_FILES_DIRECTORY_NAME] -p 12 --output-dir [SEQUENCING_FILES_DIRECTORY_NAME]/fastq_files --no-lane-splitting

######################### Sequence QC and adapter removal #########################

for sample in "${samples[@]}"
do
    :
    fastp -i "${sample}"_R1_001.fastq.gz -I "${sample}"_R2_001.fastq.gz \
    -o "${sample}"_R1_001out.fastq.gz -O "${sample}"_R2_001out.fastq.gz \
    --html "${sample}".html --json "${sample}".json

done


######################### Alignment to genome | generation of bam and BigWig files #########################

## "genome" files in .bt2 format 

for sample in "${samples[@]}"
do
    :
    echo "Starting alignment on ${sample}, go get a coffee"
    # alignment to desired genome
    bowtie2 -p 4 -x genome --very-sensitive --no-unal \
    -1 "${sample}"_R1_001out.fastq.gz -2 "${sample}"_R2_001out.fastq.gz \
    -S "${sample}"_aligned.sam
    
    # generation of bam file # filtering out bad quality reads
    samtools view -q 30 -Sb "${sample}"_aligned.sam > "${sample}"_aligned.bam
    
    samtools sort "${sample}"_aligned.bam > "${sample}"_aligned_sorted.bam
    samtools index "${sample}"_aligned_sorted.bam > "${sample}"_aligned_sorted.bam.bai
    samtools view -c "${sample}"_aligned.bam # display number of reads in bam file
    
    # generation of bigwig for visualisation pre-duplicate removal
    bamCoverage -p 4 -bs 10 --normalizeUsing BPM --exactScaling \
    -b "${sample}"_aligned_sorted.bam -o "${sample}"_aligned_sorted.bw

done


#########################  removal of mitochondrial reads ######################### 


for sample in "${samples[@]}"
do
    :
    bedtools intersect -ubam -a "${sample}"_aligned_sorted.bam -b Danio_mitoch_sorted.bam -v > "${sample}"_aligned_nucl.bam # intersect using bam files!
    samtools view -c "${sample}"_aligned_nucl.bam # print number of reads after mitochondria seq exclusion
    sambamba markdup -r -p "${sample}"_aligned_nucl.bam "${sample}"_aligned_uniq.bam # duplicate removal with Sambamba
    samtools view -c "${sample}"_aligned_uniq.bam # print number of reads after duplicate removal
    
    sambamba markdup -p "${sample}"_aligned_nucl.bam "${sample}"_ataqv.bam # for ataqv downstream processing, duplicates marked but not removed
    samtools sort "${sample}"_ataqv.bam > "${sample}"_ataqv_sorted.bam # sort and index for ataqv processing
    samtools index "${sample}"_ataqv_sorted.bam > "${sample}"_ataqv_sorted.bam.bai

done


#########################  read shifting | big wig generation ######################### 

## To account for Tn5 cutting effects
# read1: +4bp 
# read2 -5bp shift

for sample in "${samples[@]}"
do
    :
    alignmentSieve -b "${sample}"_aligned_uniq.bam --ATACshift -o "${sample}"_aligned_sieved.bam
    samtools view -c "${sample}"_aligned_sieved.bam # print number of reads after duplicate removal
    samtools sort "${sample}"_aligned_sieved.bam > "${sample}"_aligned_sieved_sorted.bam # generate indices for downstream processing with R
    samtools index "${sample}"_aligned_sieved_sorted.bam > "${sample}"_aligned_sieved_sorted.bam.bai
    
    # generate bigwig files for visualisation
    bamCoverage -p 4 -bs 10 --normalizeUsing BPM --exactScaling \
    -b "${sample}"_aligned_sieved_sorted.bam -o "${sample}"_aligned_sieved_sorted.bw 

done


#########################  calling significant peaks with macs2 ######################### 

## the following loop with produce both narrowpeaks and broadpeaks for each sample
 
for sample in "${samples[@]}"
do
    :
    mkdir "${sample}_macs2_2022_narrow"
    
    macs2 callpeak -t "${sample}"_aligned_sieved.bam -f BAM -g 1.412e9 \
    --keep-dup all --nomodel --shift -100 --extsize 200 \
    --outdir "${sample}"_macs2_2022_narrow -n "${sample}_Nar" -B -q 0.05

    mkdir "${sample}_macs2_2022_broad"

    macs2 callpeak -t "${sample}"_aligned_sieved.bam -f BAM -g 1.412e9 \
    --keep-dup all --nomodel --shift -100 --extsize 200 \
    --outdir "${sample}"_macs2_2022_broad -n "${sample}_Br" -B -q 0.05 --broad

done

#########################  TSS enrichment with ataqv ####################### 

## the dependencies here have to be reloaded because different versions are required
module add palma/2021a
module add GCC/10.3.0
module add ataqv/1.3.0


for sample in "${samples[@]}"
do
    :
    ataqv --threads 6 --peak-file "${sample}"_Br_peaks.broadPeak --name "${sample}" --ignore-read-groups \
    --metrics-file "${sample}".ataqv.json.gz --tss-file TSSRefseqZv11.bed --tss-extension 2000 \
    --autosomal-reference-file "LIST_AUTOSOMES.txt" NA "${sample}"_ataqv_sorted.bam > "${sample}".ataqv.out

    
done

# run mkarv on the JSON files to generate the interactive web viewer
mkarv [DIRECTORY_NAME] [LIST HERE ALL SAMPLES TO BE INCLUDED IN THE ANALYSIS] sample_name.ataqv.json.gz


#########################  Motif enrichment analysis with HOMER ####################### 

## Downloading homer and packages 

#perl <CURRENT_PATH>/HOMER/configureHomer.pl  -install zebrafish-p 
#perl <CURRENT_PATH>/HOMER/configureHomer.pl  -install zebrafish-o
#perl <CURRENT_PATH>/HOMER/configureHomer.pl  -install danRer11
export PATH=<CURRENT_PATH>/HOMER/bin/:${PATH}

## see DiffBind Code for generating peak/BED files (by merging replicates per condition) for use by HOMER
mkdir <OUTPUT_FOLDER/>
<CURRENT_PATH>/HOMER/bin/findMotifsGenome.pl <peak/BED file from DiffBind> danRer11 <OUTPUT_FOLDER/> -size given  -p 4




