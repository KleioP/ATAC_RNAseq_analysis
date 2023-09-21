######################### Loading dependencies and packages #########################

module add palma/2022a
module add GCC/11.3.0
module add fastp/0.23.2

module add palma/2019a
module add GCC/8.2.0-2.31.1
module add STAR-Fusion/1.8.0-Perl-5.28.1-Python-3.7.2

module add palma/2021a
module add GCC/10.3.0
module add SAMtools/1.13

module add palma/2019a
module add GCC/8.2.0-2.31.1
module add  OpenMPI/3.1.3
module add deepTools/3.3.1-Python-3.7.2
module add Subread/2.0.0
module add  MultiQC/1.9-Python-3.7.2

samples=( {SAMPLE_NAMES} )

######################### Preprocessing #########################

# where several seq runs were performed to reach 15M read samples, concatenating fastq.gz files for each sample

# loop for samples with 2 sets of reads


for sample in "${samples[@]}"
do
    :
    zcat "${sample}"_10116_4_1.fastq.gz "${sample}"_10119_4_1.fastq.gz | gzip > "${sample}"_1.fastq.gz
    zcat "${sample}"_10116_4_2.fastq.gz "${sample}"_10119_4_2.fastq.gz | gzip > "${sample}"_2.fastq.gz

done

######################### QC | duplicate removal | adapter trimming #########################

for sample in "${samples[@]}"
do
    :
    fastp -i "${sample}"_1.fastq.gz -I "${sample}"_2.fastq.gz \
    -o "${sample}"_1out.fastq.gz -O "${sample}"_2out.fastq.gz \
    --dedup --html "${sample}".html --json "${sample}".json

done


######################### Alignment to Danio rerio genome #########################

for sample in "${samples[@]}"
do
    :
    STAR --runThreadN 4 --genomeDir /home/p/petratou/myRNAseq_rep1 \
    --readFilesCommand gunzip -c --outFileNamePrefix "${sample}" \
    --readFilesIn "${sample}"_1out.fastq.gz "${sample}"_2out.fastq.gz
    echo ""${sample}" finished"
    
done

#### Genome was indexed using the following command
# STAR --runThreadN 4 --runMode genomeGenerate --genomeDir <GENOME_PATH> --genomeFastaFiles <GENOME_PATH>/GRCz11.fa


### SAM to BAM files #
### generating bigwig files for visualisation prior to read elimination

for sample in "${samples[@]}"
do
    :
    echo ""${sample}" in process, this will take a while :)"
    samtools view -q 30 -Sb "${sample}"Aligned.out.sam > "${sample}"_aligned.bam
    samtools sort "${sample}"_aligned.bam > "${sample}"_aligned_sorted.bam
    samtools index "${sample}"_aligned_sorted.bam > "${sample}"_aligned_sorted.bam.bai
    echo "Number of reads in aligned BAM file for "${sample}":"
    samtools view -c "${sample}"_aligned.bam #display number of reads in bam file. output order will be according to order of entries in the samples vector
    echo ""${sample}" finished"

done

######################### Track visualisation #########################
for sample in "${samples[@]}"
do
    :
    bamCoverage -p 4 -bs 10 --normalizeUsing BPM --exactScaling \
    -b "${sample}"_aligned_sorted.bam -o "${sample}"_RNAseq.bw

done

######################### Counting reads #########################

# Annotation file in .gtf format
#-g id ## if using .gff annotation file from ncbi

for sample in "${samples[@]}"
do
    :
    featureCounts -p -T 4 -a <ENSEMBLE ANNOTATION FILE> f -g gene_id -o "${sample}".txt "${sample}"_aligned_sorted.bam
    
done


######################### Multiqc on all log files #########################

#all subfolder under <LOG DIRECTORY> will be searched for log files
# both STAR logs and feature counts logs will be considered

multiqc <LOG DIRECTORY> -n <NAME OF OUTPUT> 
