#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter1
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-29

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/21_cardellina

# pull out header and add to filtered vcf file
grep "#" ${workdir}/03_vcf/${input_array}.g.vcf > ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf

# filter out rows that have low quality filters, genotyped sites with quality less than 20, and null alleles (* in col 4)
grep -v "#" ${workdir}/03_vcf/${input_array}.g.vcf | grep -v "LowQual" | awk '$6 >= 20 || $6 ~ /^\./' | awk '$5 !~ /*/' >> ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf

# biallelic SNPs from ingroup, for admixture, pca (no missing)
# minor allele count 2 or more, thinned to 5kbp
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --remove-indv Myioborus_miniatus_TMVB_112407 --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 2 --thin 5000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/05_structure/${input_array}_structure_5kbpthin

# biallelic SNPs from ingroup + outgroup for treemix (no missing)
# minor allele count 2 or more, thinned to 5kbp
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --max-missing 1.0 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 2 --thin 5000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/06_treemix/${input_array}_treemix_5kbpthin

# invariant and variant sites from ingroup and outgroup
# up to 4 individuals missing per snp
# for fst, dxy, heterozygosity, phylogenetics, etc.
vcftools --vcf ${workdir}/04_filtered_vcf/${input_array}.filtered.vcf --max-missing 0.8 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/07_window_stats_phylo/${input_array}


# bgzip and tabix index files that will be subdivided into windows
# directory 1
bgzip ${workdir}/07_window_stats_phylo/${input_array}.recode.vcf
#tabix
tabix -p vcf ${workdir}/07_window_stats_phylo/${input_array}.recode.vcf.gz


