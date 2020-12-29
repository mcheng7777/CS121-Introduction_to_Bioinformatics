## Usage: Rscript linreg.R

library(stringr)

args = commandArgs(trailingOnly = TRUE)
if (length(args)!=2) {
  stop("Usage: Rscript linreg.R [tped/tfam file prefix] [phenotype file]", call.=FALSE)
}

# This is the tfam file for reference
# each entry is a separate individual
# each entry contains family_id individual_id father_id mother_id sex phenotype
# NOTE: phenotype not used because we supplied an external phenotype file
t_file <- args[1]
phen_file <- args[2]
fam <- read.table(paste0(t_file,'.tfam'),header = F) 
# this is the phenotype file to regress on the snps in linear regression
# number of phenotypes = number of patients
# tfam and phentoypes contain the same ordering for patients
phenotypes <- read.table(phen_file,sep='\t', header=F)

# be sure to rearrange phenotype rows to match IDs of tfam file if necessary
phenotypes <- phenotypes[match(fam[,1],phenotypes[,1]),]

# initialize output file
# rows include chr number, SNP, p-value from linreg
outfile <- 'gwas_linreg_R.txt'
write(c('CHR','SNP','P'),file=outfile,ncolumns=3,sep='\t')


#### Linear Regression ####
# perform linear regression for each snp independantly
print('Performing linear regression for each SNP independently.')
f = file(paste0(t_file,".tped"),"r")
# initialize patient counter
p_count <- 0 
# read a line until there are no SNPs left
while (TRUE){
  # tped line: chr, SNP, genetic distance, bp, individual genotype (every 2 alleles)
  # read one line
  line = readLines(f,n=1)
  if (length(line) == 0){
    break
  }
  else {
    p_count <- p_count+1
    # split apart fields before the individual genotypes
    fields <- str_match_all(line, "(\\w*?\\d+) ")[[1]][,2]
    # extract genotypes by removing `fields` from the file line
    genotypes <- str_replace_all(str_remove(line, paste(fields,'',collapse='')), pattern = ' ', replacement = '')
    # split genotypes up into pairs of alleles (SNP)
    snps <- str_split(gsub("(.{2})", "\\1 ", genotypes),' ')[[1]]
    snps <- snps[-length(snps)]
    alleles <- unique(str_split(paste(unique(snps),collapse=''),'')[[1]])
    # encode snps as 0 1 and 2 according to genotype, there is no sense of referance allele
    encoded_snps <- str_replace(snps, pattern = paste0(alleles[1],alleles[1]), replacement = '0')
    encoded_snps <- str_replace(encoded_snps, pattern = paste0(alleles[1],alleles[2]), replacement = '1')
    encoded_snps <- str_replace(encoded_snps, pattern = paste0(alleles[2],alleles[1]), replacement = '1')
    encoded_snps <- str_replace(encoded_snps, pattern = paste0(alleles[2],alleles[2]), replacement = '2')
    # run linear regression on phenotypes and snps
    l <- summary(lm(phenotypes[,3] ~ as.numeric(encoded_snps)))
    p_val <- l$coefficients[2,4]
    # append p-value for the snp to the outfile
    write(c(fields[1],fields[2],p_val),file=outfile,ncolumns=3,sep='\t',append=T)
    # print checkpoints
    if (p_count/50000 == p_count%/%50000){
      print(paste0("Linear regression is complete for ",p_count," SNPS."))
    }
  }
   
}

close(f)


## compare p-values from GWAS and R linear regression
# R linear regression p values
print('Calculating PLINK GWAS and R p-value differences.')
r_linreg <- read.table('gwas_linreg_R.txt',header=T,sep='\t')
rpval <- r_linreg['P']
gwas <- read.table('gwas_linreg.assoc.linear',header=T)
gpval <- gwas['P']
pval_diff <- rpval-gpval
# save p-value differnce plot as png
png(filename = 'R_PLINK_pval_difference.png')
hist(x=pval_diff[,1],
     breaks=100,
     xlab = 'Difference (R p-value - PLINK p-value)',
     main="Histogram of Differences in Linear Regression P-values from R and PLINK",
     col = 'blue')
dev.off()
print('R_PLINK_pval_difference.png has been saved.')