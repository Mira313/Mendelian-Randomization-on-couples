#!/usr/bin/bash

######################
## This file take out SNPs p-values from GWAS imp
######################

######################
## Libraries
######################
print("===================")
library(tidyr)
library(data.table)
library(stringr)

######################
## Arguments from bash file
## phenotype name and localisation of GWAS p-values and sex
######################
args = commandArgs(trailingOnly=TRUE)  # for example '48_irnt' 'body' 'female'
phenotype = args[1]
folder = args[2]

input = paste0(phenotype, '.gwas.imputed_v3.both_sexes.tsv.gz')
where = paste0('/data/sgg3/data/neale_files/both_sexes/', folder, '/', input)
output = paste0(phenotype, '/', phenotype, '_pvalues.tsv')

######################
## Loading data
######################
print("1. Loading data")
data = fread(where,
        sep = "\t",
        header = TRUE,
        select = c("variant","low_confidence_variant", "pval"))

######################
## Modifying data
######################

print("2. Removing low confidence variant")
data = data[data$low_confidence_variant == FALSE,]
print(paste('nrow :', nrow(data), ', min pvalue: ', min(data$pval)))

######################
## Prunning SNPs
## Does not take X chromosome
######################
print("3. Prunning SNPs")
SNPs = data.frame()
nchr = c(1:22)
range = 5 * 10^(-8)
print(table(data$pval < range))

for (chr in nchr) {
  print(paste0('chr: ', chr))
	data_subset = subset(data, str_match(data$variant, "\\s*(.*?)\\s*:")[,2] == chr)    # taking chromosome of SNP
	pos = as.vector(as.numeric(str_match(data_subset$variant, ":\\s*(.*?)\\s*:")[,2]))  # taking position of SNP
	pval_table = data_subset$pval < range                                                      
	while (TRUE %in% pval_table) {
		n = which(data_subset$pval == min(data_subset$pval))                              # which row has the lowest pval
		SNPs = rbind(SNPs, data_subset[n,])                                               # adding SNP to list
		data_subset = data_subset[!(pos <= pos[n]+500000 & pos >= pos[n]-500000),]        # removing SNP and neighbours from data
		pos = pos[!(pos <= pos[n]+500000 & pos >= pos[n]-500000)]
		pval_table = data_subset$pval < range                                             # update condition
	}
}

######################
## Writing file
######################
print("4. Writing in file")
print(paste0('nrow: ', nrow(SNPs)))
write.table(SNPs, output, sep ="\t", row.names = FALSE) 

######################
## END
######################
print("===================")
