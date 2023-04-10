# start interactive session
interactive -p quanah

# move to directory
workdir=/lustre/scratch/jmanthey/21_cardellina
cd ${workdir}/05_structure

# no Z chromosome
rm CM027535.1_structure_5kbpthin.recode.vcf

# combine all chromosome vcfs into one vcf for the two different datasets
grep "#" CM027537.1_structure_5kbpthin.recode.vcf > 5kbpthin.vcf
for i in $( ls *structure_5kbpthin* ); do grep -v "#" $i >> 5kbpthin.vcf; done
rm *recode*

# make chromosome map for the vcf
grep -v "#" 5kbpthin.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the combined vcf
vcftools --vcf 5kbpthin.vcf  --plink --chrom-map chrom_map.txt --out 5kbpthin 

# convert  with plink
plink --file ${workdir}/05_structure/5kbpthin --recode12 --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink

# run pca
plink --file ${workdir}/05_structure/5kbpthin_plink --pca --allow-extra-chr --out ${workdir}/05_structure/5kbpthin_plink_pca

# run admixture
for K in 1 2 3 4 5; do admixture --cv ${workdir}/05_structure/5kbpthin_plink.ped $K; done
