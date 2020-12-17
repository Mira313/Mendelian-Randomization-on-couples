#!/usr/bin/bash

######################
## This file does the correlation between the 20 first principal components 
## and the phenotypes of both sexes combined and the correlation between sexes for each PC
######################

######################
## Libraries
######################
print("===================")
library(ggplot2)
library(gplots)
library(ggExtra)

######################
## Performing correlation between PCs and phenotypes
######################
print("===================")
PCs = read.table('PCs.txt', sep = ' ', header = T)

female.Id = read.table('femaleIds.txt', sep = ' ', header = T)
male.Id = read.table('maleIds.txt', sep = ' ', header = T)
colnames(female.Id) = c('Id', 'age')
colnames(male.Id) = colnames(female.Id)

phenotypes = read.table('phenotypes.txt', sep = ' ', header = T)
phenotypes$Neale_pheno_ID = as.character(phenotypes$Neale_pheno_ID)

PCs$Id = NULL

PCs.female = PCs[1:(nrow(PCs)/2),]
PCs.male = PCs[((nrow(PCs)/2)+1):nrow(PCs), ]

cor.pheno_PC = as.data.frame(matrix(NA, nrow = nrow(phenotypes), ncol = 20))
colnames(cor.pheno_PC) = colnames(PCs)
rownames(cor.pheno_PC) = phenotypes$Neale_pheno_ID

for (pheno in rownames(pheno)) {
  loc = paste0(pheno, '/', pheno)
  print(loc)
  file = paste0(loc, '_female_pheno.txt')
  female.pheno = read.table(file, sep = ' ', header = T)
  file = paste0(loc, '_male_pheno.txt')
  male.pheno = read.table(file, sep = ' ', header = T)
  
  pheno.data = rbind(female.pheno, male.pheno)
  
  for (pc in colnames(pheno)) {
    cor.pheno_PC[pheno,pc] = cor.test(PCs[,pc], pheno.data[,1])$estimate
  }
}
rownames(cor.pheno_PC) = phenotypes$description
head(cor.pheno_PC)
write.table(cor.pheno_PC, 'PC_pheno.txt', col.names = T, row.names = T)

######################
## Doing heatmap
######################
cor.pheno_PC = as.matrix(cor.pheno_PC)

pdf('../plots/heatmap.pdf', paper = "USr", width = 11, height = 8)
heatmap.2(cor.pheno_PC, Colv = F, dendrogram = 'row', scale = 'none', trace = 'none',
          density.info = 'histogram', key = T, key.xlab = 'Correlation', keysize = 0.9,
          key.title = '', margins = c(5,15.5), cexRow = 0.5 + 1/log10(nrow(r)),
          cexCol = 0.4 + 1/log10(nrow(r)), rowsep = c(10, 23))
dev.off()

######################
## Performing correlation between sexes for each PC
######################
cor.sex_PC = as.data.frame(matrix(NA, nrow = ncol(PCs), ncol = 3))
colnames(cor.sex_PC) = c('cor', 'conf.inf', 'conf.sup')
rownames(cor.sex_PC) = colnames(PCs)

for (pc in colnames(PCs)) {
  r = cor.test(PCs.female[ ,pc], PCs.male[, pc])
  cor.sex_PC[pc,] = c(r$estimate, r$conf.int[1], r$conf.int[2])
}
write.table(cor.sex_PC, 'PC_sex.txt', sep = ' ', row.names = T, col.names = T)

######################
## Doing barplot
######################
results$Id = rownames(results)
results$Id = factor(results$Id, levels = results$Id)

pdf('../plots/barplot.pdf')
ggplot(data=results, aes(x=Id, y=cor)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=conf.inf, ymax=conf.sup), width=.2,
                position=position_dodge(.9)) +
  theme_classic(base_size = 9) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(position = "right") +
  ylab('Correlation \n between sexes') +
  removeGridX()
dev.off()

######################
## END
######################
print("===================")

