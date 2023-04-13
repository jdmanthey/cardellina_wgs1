interactive -p quanah

workdir=/lustre/scratch/jmanthey/21_cardellina/06_treemix

cd ${workdir}

# move all the files into a single vcf
grep "#" CM027507.1_treemix_5kbpthin.recode.vcf > treemix.all.vcf
for i in $( ls *recode.vcf ); do
	grep -v "#" $i >> treemix.all.vcf;
done

# bgzip and index the vcf files
bgzip -c ${workdir}/treemix.all.vcf > ${workdir}/treemix.all.vcf.gz

tabix -p vcf ${workdir}/treemix.all.vcf.gz

# run bcftools to simplify the vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/treemix.all.vcf.gz > ${workdir}/treemix.simple.vcf

# convert in R for treemix
# load R
module load intel R
# start R
R

# vcf to treemix function


vcf_to_treemix <- function(input, output, popmap) {
	options(scipen=999)
	
	# read in input file
	x <- read.table(input)
	
	# read in popmap 
	popmap2 <- read.table(popmap, stringsAsFactors=F, header=T)
	populations <- unique(popmap2[,2])
	
	# replace all phased genotypes with unphased genotypes
	for(a in 5:ncol(x)) {
		x[,a] <- gsub("\\|", "/", x[,a])
	}
	
	# check each snp to make sure it is in each population
	keep <- rep(TRUE, nrow(x))
	for(a in 1:nrow(x)) {
		a_rep <- as.character(x[a,5:ncol(x)])
		for(b in 1:length(populations)) {
			if(keep[a] == TRUE) { # only check in this pop if the snp still included
				b_rep <- unique(a_rep[popmap2[,2] == populations[b]])
				b_rep <- b_rep[b_rep != "./."]
				if(length(b_rep) == 0) {
					keep[a] <- FALSE
				}
			}			 
		}
	}
	x <- x[keep,]
	
	# write output per snp
	write(populations, file=output, sep=" ", ncolumns=length(populations), append=F)
	for(a in 1:nrow(x)) {
		a_rep <- as.character(x[a,5:ncol(x)])
		a_output <- c()
		for(b in 1:length(populations)) {
			b_rep <- a_rep[popmap2[,2] == populations[b]]
			b_allele1 <- length(b_rep[b_rep == "0/0"]) * 2 + length(b_rep[b_rep == "0/1"])
			b_allele2 <- length(b_rep[b_rep == "1/1"]) * 2 + length(b_rep[b_rep == "0/1"])
			a_output <- c(a_output, paste0(b_allele1, ",", b_allele2))
		}
		write(a_output, file=output, sep=" ", ncolumns=length(populations), append=T)
	}	
}

# convert to treemix format
vcf_to_treemix("treemix.simple.vcf", "cardellina.treemix", "popmap_treemix.txt")

q()

gzip cardellina.treemix


## run treemix 


# run with up to three migration edges:
src/treemix -i input/cardellina.treemix.gz -root outgroup -o output/cardellina.treemix
src/treemix -i input/cardellina.treemix.gz -m 1 -g output/cardellina.treemix.vertices.gz output/cardellina.treemix.edges.gz -o output/cardellina_m1.treemix
src/treemix -i input/cardellina.treemix.gz -m 1 -g output/cardellina_m1.treemix.vertices.gz output/cardellina_m1.treemix.edges.gz -o output/cardellina_m2.treemix
src/treemix -i input/cardellina.treemix.gz -m 1 -g output/cardellina_m2.treemix.vertices.gz output/cardellina_m2.treemix.edges.gz -o output/cardellina_m3.treemix


#bootstraps over 200s snps without migration edges
for i in {1..100}; do
    src/treemix -i input/cardellina.treemix.gz -root outgroup -bootstrap -k 200 -o output/$i.treemix
done;

# unzip the tree files
for i in { ls *treeout.gz }; do
    gzip -d $i
done;

# in R:
#summarize bootstraps
x <- list.files(pattern="*treeout")
for(a in 1:length(x)) {
if (a==1) {
output <- scan(x[a], what="character")[1]
} else {
output <- c(output, scan(x[a], what="character")[1])
}
}
write(output, file="cardellina_treemix_bootstraps.trees", ncolumns=1)

# in bash
# summarize bootstraps
sumtrees.py --output=cardellina_treemix_summed.tre --min-clade-freq=0.01 cardellina_treemix_bootstraps.trees




