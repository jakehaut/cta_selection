# cta_selection 
usage: ./test_bash.sh [snp list file, space deliminated] [plink file name] [output file name] [number of repitions] [epsilon, similar snp frequency tolerance] [bounds to use - sigma, or 95%CI otherwise]

example: ./cta_selection_script.sh height_snps_top_nomissing.txt allchrmsmerged_eur50subs_nodups heightsnpstop_eur50_100r0.01e 100 0.01 ci

You need a file system set up as follows:
-analysis folder, containing: cta_selection_script.sh, plot_relmat.py, generate_cta_snps.py and the snp list file
-data folder, containing: .bim/.bed/.fam plink files for the whole genome, subsetted to contain the individuals you want to look at
-tools folder, at your home directory (or wherever '~/' points to) containing plink v1.9
