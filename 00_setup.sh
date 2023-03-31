# define main working directory
workdir=/lustre/scratch/jmanthey/21_cardellina

# make output directories
cd ${workdir}

mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_filtered_vcf
mkdir 05_structure
mkdir 06_treemix
mkdir 07_window_stats_phylo
mkdir 08_sfs
mkdir 10_filter
mkdir 12_filter_vcf

# index reference
interactive -p quanah

module load intel java bwa

samtools faidx GCA_001746935.2_mywa_2.1_genomic.fna
bwa index GCA_001746935.2_mywa_2.1_genomic.fna
java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/GCA_001746935.2_mywa_2.1_genomic.fna O=/home/jmanthey/references/GCA_001746935.2_mywa_2.1_genomic.dict
