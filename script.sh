
### This is the main script to run all commands of the project
# Usage: bash script.sh [bed/bim/fam/tped/tfam prefix] [phenotype filename] [UK Biobank GWAS .tsv file]
# check appropriate number of arguments
if [ "$#" -ne 3 ];
then
    echo "Usage: bash script.sh [bed/bim/fam file prefix] [phenotype filename] [UK Biobank GWAS .tsv file]"
    exit 1
fi

snp_data_file_prefix=$1
phenotype_file=$2
UK_Biobank_HDL_file=$3

# check if files exist in directory
if [ ! -f ${snp_data_file_prefix}.bed ];
then
    echo "${snp_data_file_prefix}.bed does not exist"
    exit 1
fi
if [ ! -f ${snp_data_file_prefix}.bim ];
then
    echo "${snp_data_file_prefix}.bim does not exist"
    exit 1
fi
if [ ! -f ${snp_data_file_prefix}.fam ];
then
    echo "${snp_data_file_prefix}.fam does not exist"
    exit 1
fi
if [ ! -f ${phenotype_file} ];
then
    echo "${phenotype_file} does not exist"
    exit 1
fi
if [ ! -f ${UK_Biobank_HDL_file} ];
then
    echo "${UK_Biobank_HDL_file} does not exit"
    exit 1
else
    if [[ ${UK_Biobank_HDL_file} != *'.tsv' ]];
    then
        echo "${UK_Biobank_HDL_file} must be a tsv file"
        exit 1
    fi
fi

# Step 1: PLINK commands - GWAS and file conversion
echo "Now running PLINK Commands: GWAS and file conversion\n\n"
bash gwas.sh ${snp_data_file_prefix} ${phenotype_file}

# check if Step 1 output files are written
if [ ! -f ${snp_data_file_prefix}.tped ];
then
    echo "${snp_data_file_prefix}.tped file was not written"
    exit 1
fi
if [ ! -f ${snp_data_file_prefix}.tfam ];
then
    echo "${snp_data_file_prefix}.tfam file was not written"
    exit 1
fi
if [ ! -f gwas_linreg.assoc.linear ];
then
    echo "gwas_linreg.assoc.linear file was not written"
    exit 1
fi

# Step 2: Linear regression in R using the tped and tfam files
# check if necessary files exist

echo "Now running linear regression with R with linreg.R. Outputs are:
- gwas_linreg_R.txt - chromosome, rsid, and p-value of the snp from linear regression
- R_PLINK_pval_difference.png - histogram of p-value difference between R linear regression and PLINK GWAS for each SNP \n\n"

Rscript linreg.R ${snp_data_file_prefix} ${phenotype_file}

# Step 3: Manhattan Plot
echo "Now creating manhattan plot with qqman.r. Outputs are:
- manhattan_plot.png - a manhattan plot that plots the -log10 pvalue for each SNP\n\n"
Rscript qqman.r

# Step 4: Compare minimal p-value SNPs of each chromosome to UK Biobank HDL Choleterol GWAS
echo "Now comparing minimal p-value SNPs of each chromosome from PLINK GWAS to UK Biobank HDL Cholesterol GWAS. Outputs are:
- PLINK_HDL_SNP.txt - table of minimal p-value snps of each chromosome from PLINK GWAS and the p-values from PLINK and UKBiobank GWAS\n\n"
Rscript HDL_analysis.R ${UK_Biobank_HDL_file}


