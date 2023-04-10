# Part 1

#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=cat_vcf
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

chr_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_cat.txt | tail -n1 )

chr_start=${chr_array}a.g.vcf

chr_output=${chr_array%__}.g.vcf

grep "#" $chr_start > $chr_output

for i in $( ls $chr_array*vcf ); do grep -v "#" $i >> $chr_output; done


# Part 2

# remove all partial chromosome vcfs
rm *__*

# combine those chromosomes that were split
grep -v "#" CM027507.1b.g.vcf >> CM027507.1.g.vcf
grep -v "#" CM027508.1b.g.vcf >> CM027508.1.g.vcf
grep -v "#" CM027509.1b.g.vcf >> CM027509.1.g.vcf

# remove b files
rm CM027507.1b.g.vcf
rm CM027508.1b.g.vcf
rm CM027509.1b.g.vcf





