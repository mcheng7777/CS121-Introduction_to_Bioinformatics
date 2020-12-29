
# Usage: bash gwas.sh [bed/bim/fam file prefix] [phenotype filename]
# run PLINK commands
# 1) gwas on phenotypes to obtain SNP marginal effects and p-values
# 2) recode the bed files to get tped and tfam files

# check appropriate number of arguments
if [ "$#" -ne 2 ];
then
    echo "Usage: gwas.sh [bed/bim/fam file prefix] [phenotype filename]. Requires 2 arguments"
    exit 1
fi


# set bed files and phenotype files directories
bin_files=$1
phen_file=$2

echo "Running GWAS. Outputs are:
- gwas_linreg.assoc.linear - contains summary statistics for each SNP, including the effect size and p-value"

# run plink linear regression (GWAS)
plink \
    --bfile $bin_files \
    --linear \
    --pheno $phen_file \
    --out "gwas_linreg"

echo "Finished GWAS"

# convert bed files to transposed tped/tfam files
echo "Converting bed files to tped and tfam files. Outputs are:
- snpdata.tped - transposed ped file, contains SNP id and genotype for each individual
- snpdata.tfam - contains ID information for each individual"
plink \
    --bfile $bin_files \
    --recode transpose \
    --out $bin_files
echo "Finished converting"
