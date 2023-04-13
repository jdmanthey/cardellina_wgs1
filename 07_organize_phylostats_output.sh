# headers
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_theta.txt
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_pi.txt
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_tajima.txt
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_dxy.txt
grep 'pop1' CM027507.1__100000001__100050000__stats.txt > ../window_fst.txt

# add the relevant stats to each file
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *txt ); do grep 'theta' $i >> ../window_theta.txt; done
for i in $( ls *txt ); do grep 'pi' $i >> ../window_pi.txt; done
for i in $( ls *txt ); do grep 'Tajima_D' $i >> ../window_tajima.txt; done
for i in $( ls *txt ); do grep 'Dxy' $i >> ../window_dxy.txt; done
for i in $( ls *txt ); do grep 'Fst' $i >> ../window_fst.txt; done

# combine all trees
cat *bipart* >> ../50kbp.trees
