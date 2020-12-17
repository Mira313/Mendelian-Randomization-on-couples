#!/usr/bin/bash

######################
## This file does prunning on SNPs pvalues from GWAS
## And take out SNPs for males and females for each phenotype
######################

# 'body' ('48_irnt' '49_irnt' '50_irnt' '21001_irnt' '21002_irnt' '23099_irnt' '23116_irnt' '23117_irnt' '23124_irnt' '23125_irnt' '23128_irnt' '23129_irnt' '23105_irnt')
# 'lifestyle' ('1558' '845' '20016_irnt' '3456' '874_irnt' '924')
# 'diet' ('1289' '1299' '1309' '1319' '1329' '1339' '1349' '1359' '1369' '1379' '1389' '1408' '1438_irnt' '1458' '1478' '1488_irnt' '1498' '1528' '1548')

######################
## Libraries
######################
print("===================")
library(rbgen)
library(data.table)
library(tidyr)
library(stringr)

print(format(Sys.time(), "%X"))

######################
## Prunning
######################
print('----------PRUNNING STEP-----------')

######################
## Arguments from bash file
## phenotype name and localisation of GWAS p-values
######################
args = commandArgs(trailingOnly=TRUE)  # for example '48_irnt' 'body'
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
## Taking out SNPs
######################
print('----------TAKING OUT SNPs STEP-----------')

######################
## Loading Datas
######################

print("1. Reading Ids")
fID = read.table("femaleIds.txt", sep = " ", header = TRUE)       # taking IDs of females
mID = read.table("maleIds.txt", sep = " ", header = TRUE)         # taking IDs of males
d = read.table("/data/sgg3/data/UKBB/imp/ukb1638_imp_chr1_v2_s487398.sample", sep = " ", header = TRUE) 
d = d[2:nrow(d),]

print("3. Comparing Ids")
whichFId = (d$ID_1 %in% fID$femaleId)
whichMId = (d$ID_1 %in% mID$maleId)

print("4. Reading rsids")
rsids = fread("variants.tsv",
              sep = "\t",
              header = TRUE,
              select = c('variant', 'rsid'))

######################
## Reading rsids
######################

nchr = c(1:22)
sex = c('male', 'female')

print("5. Reading pvalues")

print("Comparing rsids")
whichSid = (rsids$variant %in% SNPs$variant)
my_rsids = rsids$rsid[whichSid]
SNPs$rsid = my_rsids
print(paste0('nrow pval: ', nrow(SNPs)))

for (s in sex) {
	print(paste0('-------------- ', s))
  
	print("Reading bgen file")
  first = TRUE
	if (s == 'male') {
		whichId = whichMId
	} else if (s == 'female') {
		whichId = whichFId
	}
	DATA = data.frame(matrix(NA, nrow = nrow(fID), ncol = 1))
  		
	n = 0
	total_rsids = c()
	for (chr in nchr) {
	  my_rsids2 = my_rsids[str_match(SNPs$variant, "\\s*(.*?)\\s*:")[,2] == chr]
		total_rsids = append(total_rsids, my_rsids2)
    
		file = paste0('/data/sgg3/data/UKBB/imp/_001_ukb_imp_chr', chr, '_v2.bgen')
    bgen = bgen.load(file, rsids = my_rsids2)
    data = bgen$data[, whichId, 2] + 2 * bgen$data[, whichId, 3]
      	
		if (is.null(dim(data)) == TRUE) {
			DATA = cbind(DATA, data)
			names(DATA)[length(names(DATA))] = my_rsids2[1]
		} else {
			DATA = cbind(DATA, t(data))
		}
  }
  
  print("Matching IDs")
	DATA[1] = NULL
	output = paste(phenotype, '/', phenotype, '_', s, '_geno.tsv', sep = '')
	if (s == 'male') {
    D = data.frame(ID = d$ID_1[whichMId], DATA)
    final_data = D[match(mID$maleId , D$ID), ]
  } else if (s == 'female') {
    D = data.frame(ID = d$ID_1[whichFId], DATA)
    final_data = D[match(fID$femaleId , D$ID), ]
  }
	print(paste0('Data dimension: ', dim(final_data)))
	write.table(final_data, output, sep = '\t', row.names=FALSE)	
}

######################
## END
######################
print("===================")

