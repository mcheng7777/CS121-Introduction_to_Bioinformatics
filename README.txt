CS 121 Project #3

In order to run the entire pipeline at once, run script.sh, which will call the commands for all the smaller scripts to complete the entire project. Here is how you use it:

bash script.sh [bed/bim/fam/tped/tfam prefix] [phenotype filename] [UK Biobank GWAS .tsv file] 

Example: if you have in your directory
- snpdata.bim, snpdata.bed, snpdata.fam, phenotype.txt, and HDL.tsv
Run the following code:

bash script.sh snpdata phenotype.txt HDL.tsv