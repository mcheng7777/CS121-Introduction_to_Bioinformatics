
library(dplyr)
library(stringr)


args = commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: Rscript HDL_analysis.R [UK Biobank gwas file]", call.=FALSE)
}

HDL_file <- args[1]
# read in PLINK GWAS and HDL GWAS file
print('Reading in PLINK GWAS and HDL UK Biobank GWAS')
gwas <- read.table('gwas_linreg.assoc.linear',header=T)
HDL <- read.table(HDL_file,header=T)
HDL <- HDL[,c('variant','pval')]
HDL$variant <- str_match(HDL$variant,'\\d+:\\d+')

# selecting minimal p-value SNPs from each chromosome of PLINK GWAS
print("Selecting minimal p-value SNPs.")
gwas_chrom <- as.data.frame(gwas %>% group_by(CHR) %>% slice(which.min(P)))

# combine chr and bp of gwas
g <- paste(gwas_chrom$'CHR',gwas_chrom$'BP',sep=':')

# extract corresponding SNPs from HDL GWAS
HDL_chrom <- HDL[which(HDL$'variant' %in% g),]
HDL_chrom$'order' <- match(HDL_chrom$'variant',g)
HDL_chrom <- HDL_chrom %>% arrange(order)

# create dataframe with both PLINK GWAS and UK Biobank GWAS p-values
print('Writing p-value dataframe.')
gwas_chrom$'HDL_P' <- HDL_chrom$'pval'
colnames(gwas_chrom)[9] <- 'PLINK_P'
write.table(gwas_chrom[,c('SNP','PLINK_P','HDL_P')],file = 'PLINK_HDL_SNP.txt',sep='\t',col.names = T,row.names=F,quote=F)
print('PLINK_HDL_SNP.txt has been saved.')
