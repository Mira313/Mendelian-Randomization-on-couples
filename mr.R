#!/usr/bin/bash

######################
## This file take out SNPs p-values from GWAS imp
######################

######################
## Libraries
######################
print("===================")
library(tidyverse)
library(data.table)
library(devtools)
library(TwoSampleMR)
library(stringr)

rsids = fread("variants.tsv", sep = "\t", header = TRUE, select = c('variant', 'rsid', 'AF', 'alt'))

body = c('49_irnt', '50_irnt', '21001_irnt', '21002_irnt', '23099_irnt', '23116_irnt', '23117_irnt', '23124_irnt', '23125_irnt', '23128_irnt', '23129_irnt', '23105_irnt')
lifestyle = c('1558', '845', '20016_irnt', '3456', '874_irnt', '924')
diet = c('1289', '1299', '1309', '1319', '1329', '1339', '1349', '1359', '1369', '1379', '1389', '1408', '1438_irnt', '1458', '1478', '1488_irnt', '1498', '1528', '1548')

args = commandArgs(trailingOnly=TRUE)  # for example '48_irnt' 'body'
loc = args[1]
if (loc == 'body') {
  phenotypes = body
} else if (loc == 'lifestyle') {
  phenotypes = lifestyle
} else if (loc == 'diet') {
  phenotypes = diet
}

##################
## Performing analysis
######################
sex1 = c('male', 'female') # sex for genotype
sex2 = c('female', 'male') # sex for the phenotype
for (sex in 1:2) {
  s1 = sex1[sex]
  s2 = sex2[sex]

  for (p in phenotypes) {
    print(paste0('=============', p, '============'))
    pheno = paste0(p, '/', p)
    # loading which SNPs
    file = paste0(pheno, '_pvalues.tsv')
    pval = fread(file, sep = "\t", header = TRUE, select = 'variant')
    # loading genotype file
    file = paste0(pheno, '_', s1, '_geno.tsv')
    geno = read.table(file, sep = '\t', header = TRUE)
    # loading phenotype file
    file = paste0(pheno, '_', s2, '_pheno.txt')
    pheno = read.table(file, sep = ' ', header = TRUE)
    # loading beta and se estimate for exposure
    file = paste0('/data/sgg3/data/neale_files/', s1, '/', loc, '/', p, '.gwas.imputed_v3.', s1, '.tsv.gz')
    data = fread(file, sep = "\t", header = TRUE, select = c("variant","beta", "se"))
    # loading Id and age of sex 1
    file = paste0(s1, 'Ids.txt')
    Id1 = read.table(file, sep = ' ', header = T)
    # loading Id and age of sex 2
    file = paste0(s2, 'Ids.txt')
    Id2 = read.table(file, sep = ' ', header = T)
    # loading PCs
    file = '/data/sgg3/data/UKBB/geno/ukb_sqc_v2.txt'
    genetic_pcs = fread(file, sep = ' ', select = c(26:45)) # read only the first 20 PCs
    colnames(genetic_pcs) = paste0(rep('PC', 20), 1:20)
    
    age = Id2[,2]
    
    rsids2 = rsids[rsids$variant %in% pval$variant,]      # take rsids that were chosen by their GWAS p-values
    rsids2 = rsids2[match(pval$variant, rsids2$variant),] # match these variants
    
    df = data[data$variant %in% rsids2$variant,]          # retrieve variants of chosen SNPs
    df = df[match(rsids2$variant, df$variant),]           # match these variants
    df$SNP = rsids2$rsid                                  # add rsids
    
    # filtering PCs with IDs
    d = read.table("/data/sgg3/data/UKBB/plink/_001_ukb_cal_chr1_v2.fam", sep = " ", header = FALSE) 
    whichId = (d[,1] %in% Id1[,1])
    genetic_Id = d[whichId, 1]
    genetic_pcs = genetic_pcs[whichId,]
    genetic_pcs$Id = genetic_Id
    genetic_pcs = genetic_pcs[match(Id1[,1], genetic_pcs$Id),]
    geno$ID = NULL # remove ID column of genotype data
    
    print(paste0('nrow df: ', nrow(df), ' ncol geno: ', ncol(geno), ' nrow PCs: ', nrow(genetic_pcs)))
    
    ######################
    ## Removing empty columns
    ######################
    tonull_i = c()
    tonull_n = c()
    for (i in 1:ncol(geno)) {
      if (var(geno[i]) == 0) {
        tonull_i = append(tonull_i, i)
        tonull_n = append(tonull_n, colnames(geno)[i])
        
      }
    }
    print(paste0('low variance SNP: ', length(tonull_i)))
    if (length(tonull_i) > 0) {
      for (i in 1:length(tonull_n)) {
        if (substring(tonull_n[i],1,1)=="X") {
          tonull_n[i] = substr(tonull_n[i],2,nchar(tonull_n[i]))
          tonull_n[i] = gsub("\\.", ":", tonull_n[i])
        }
        df = df[df$SNP != tonull_n[i],]
        rsids2 = rsids2[rsids2$rsid != tonull_n[i],]
      }
      geno[tonull_i] = NULL
    }
    
    print(paste0('nrow df: ', nrow(df), ' ncol geno: ', ncol(geno)))
    
    ######################
    ## checking if rs are in double
    ## removing the one with the lowest variance
    ######################
    det = colnames(geno[grep('\\.', colnames(geno))])
    det = det[grep('rs', det)]
    print(paste0('Number of SNP in double: ', length(det)))
    for (snp in det) {
      occurence = which(colnames(geno) == snp)
      origin = occurence - 1
      print(paste(colnames(geno)[occurence], colnames(geno)[origin], sep = ' - '))
      print(paste(var(geno[,occurence]), var(geno[,origin]), sep = ' - '))
      if (var(geno[,occurence]) > var(geno[,origin])) {
        geno[,origin] = geno[,occurence]
        
      }
      geno[,occurence] = NULL
      print(var(geno[,origin]))
    }
    
    print(paste0('nrow df: ', nrow(df), ' ncol geno: ', ncol(geno)))
    nrsid = ncol(geno)
    
    ######################
    ## Linear regression
    ######################
    geno = cbind(geno, age, genetic_pcs)
    model = summary( lm( pheno[ , 1 ] ~ . + I(age^2), data = geno ) )
    
    ######################
    ## Creating outcome and exposure datas
    ######################
    outcome = data.frame('SNP' = rownames(model$coefficients)[2:(nrsid+1)], model$coefficients[2:(nrsid+1),1:2], row.names = NULL)
    colnames(outcome) = c('SNP', 'beta', 'se')
    outcome$SNP = as.character(outcome$SNP)
    for (i in 1:nrow(outcome)) {
      if (substring(outcome$SNP[i],1,1)=="X") {
        outcome$SNP[i] = substr(outcome$SNP[i],2,nchar(outcome$SNP[i]))
        outcome$SNP[i] = gsub("\\.", ":", outcome$SNP[i])
      }
    }
    
    outcome = outcome[match(df$SNP, outcome$SNP),]
    
    exposure = data.frame(df$SNP, df$beta, df$se, rsids2$AF)
    exposure$effect_allele = as.character(rsids2$alt)
    colnames(exposure) = c('SNP', 'beta', 'se', 'eaf', 'effect_allele')
    outcome$eaf = exposure$eaf
    outcome$effect_allele = exposure$effect_allele
    
    ######################
    ## Format, harmonise and mr
    ######################
    exp_data = format_data(exposure, type = 'exposure')
    out_data = format_data(outcome, type = 'outcome')
    har_data = harmonise_data(exp_data, out_data, action = 1)
    het_data = mr_heterogeneity(har_data)
    het_data$outcome = paste0(p, '_', s2)
    het_data$exposure = paste0(p, '_', s1)
    
    result = mr(har_data, method_list=c("mr_ivw"))
    result$outcome = paste0(p, '_', s2)
    result$exposure = paste0(p, '_', s1)
    
    print(result)
    # mr_scatter_plot(result, har_data)
    
    ######################
    ## Writing file
    ######################
    write.table(result, 'results_final.tsv', sep = '\t', row.names=F, col.names = F, append = T)
    write.table(het_data, 'heterogeneity_final.tsv', sep = '\t', row.names=F, col.names = F, append = T)
  }
}
######################
## END
######################
print("===================")


