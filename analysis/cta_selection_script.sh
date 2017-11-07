#!/bin/bash

#usage: ./test_bash.sh [snp list file, space deliminated] [plink file name] [output file name] [number of repitions] [epsilon, similar snp frequency tolerance] [bounds to use - sigma, or 95%CI otherwise]
#example: ./cta_selection_script.sh eur50_testsnps.txt chrm10_nodups_50eur chrm10_50eur_100r0.005e_testing 100 0.005 sigma
#You need a file system set up as follows:
#analysis folder, containing this script, plot_relmat.py, generate_cta_snps.py and the snp list
#data folder, with the .bim/.bed/.fam plink files for the whole genome, for the individuals you want to look at
#tools folder, at your home directory (or wherever '~/' points to) containing plink v1.9
start_time=`date +%s`
mkdir temp
echo "made temp folder"
cd temp/

echo $2
OBSSNPS=$(cat ../$1)
~/tools/./plink --bfile ../../data/$2 --freq --snps $OBSSNPS --out $2_interest > ../$2_log.txt

~/tools/./plink --bfile ../../data/$2 --freq --out $2 > ../$2_log.txt

echo "made freq file"

python ../generate_cta_snps.py ../$1 $2 $4 $5

echo "made replicate snp files"

SNPS=$(cat refsnps_observed.txt)
~/tools/./plink --bfile ../../data/$2 --snps $SNPS --make-rel square --out $3_interest >> ../$2_log.txt

for (( i = 0; i < $4; i++ )); do
	SNPS=$(cat refsnps$i.txt)
	~/tools/./plink --bfile ../../data/$2 --snps $SNPS --make-rel square --out $3_$i >> ../$2_log.txt
done
echo "made relationship matrices"

cd ..

python plot_relmat.py temp/$3 $4 $3 $6

rm -r temp/

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

