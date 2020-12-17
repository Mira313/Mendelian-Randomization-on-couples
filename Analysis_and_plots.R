#!/usr/bin/bash

######################
## This file analyse differences 
# between causation and correlation 
# between male-to-female and female-to-male causal effects
# does plots of these
######################

######################
## Libraries
######################
print("===================")
library(meta)
library(ggplot2)
library(ggrepel)

##################
## Performing analysis
######################
ivw = read.table('../data_files/ivw.txt', sep = '\t', header = T) # results from the mr

pval_difference  =  2 * pnorm( -abs( (ivw$b_m_f - ivw$b_f_m) / sqrt(ivw$se_m_f^2 + ivw$se_f_m^2) ))
p = p.adjust(pval_difference, method = 'bonferroni')
ivw$t.test = pval_difference
ivw$t.test.adjust = p

meta.data = as.data.frame(matrix(NA, nrow = nrow(ivw), ncol = 2))
colnames(meta.data) = c('TE.fixed', 'seTE.fixed')

for (i in 1:nrow(ivw)) {
  r = metagen(c(ivw$b_m_f[i], ivw$b_f_m[i]), c(ivw$se_m_f[i], ivw$se_f_m[i]))
  meta.data$TE.fixed[i] = r$TE.fixed
  meta.data$seTE.fixed[i] = r$seTE.fixed
}

# meta.data$description = ivw$description
# write.table(meta.data, '../data_files/meta_data.txt', sep = ' ', col.names = T, row.names = F)
# write.table(ivw, '../data_files/ivw.txt', sep = '\t', col.names = T, row.names = F)

##################
## Doing plots
# causation on correlation between sexes for each phenotype
# male-to-female on female-to-male causal effects
######################
meta.data = cbind(meta.data, ivw[, 1:4])
pdf('../plots/causation_correlation.pdf', paper = "USr", width = 11, height = 8)
ggplot(data = meta.data, aes(x = correlations, y = TE.fixed)) +
  labs(x = 'Phenotype correlation between sexes',
       y = 'Meta-analysed causal effect',
       color = 'Genotype sex') +
  geom_errorbar(data = meta.data, aes(ymin=TE.fixed-(1.96*seTE.fixed), ymax=TE.fixed+(1.96*seTE.fixed)), width=0.01, color = 'gray55') +
  geom_errorbarh(data = meta.data, aes(xmin = conf.inf, xmax = conf.sup), height = 0.01, color = 'gray55') +
  geom_abline(intercept = 0, slope = 1, color = 'red2') +
  geom_point() +
  geom_text_repel(size = 4, label = meta.data$desc, show.legend = FALSE) +
  theme_grey(base_size = 14)
dev.off()

pdf('../plots/male_female.pdf', paper = "USr", width = 11, height = 8)
ggplot(data = ivw, aes(x = b_f_m, y = b_m_f)) +
  labs(x = 'Female-to-male causal effect',
       y = 'Male-to-female causal effect') +
  geom_errorbar(data = ivw, aes(ymin=b_m_f-(1.96*se_m_f), ymax=b_m_f+(1.96*se_m_f)), width=0.01, color = 'gray55') +
  geom_errorbarh(data = ivw, aes(xmin=b_f_m-(1.96*se_f_m), xmax=b_f_m+(1.96*se_f_m)), height = 0.01, color = 'gray55') +
  geom_abline(intercept = 0, slope = 1, color = 'red2') +
  geom_point() +
  geom_text_repel(size = 4, label = ivw$description, show.legend = FALSE) +
  theme_grey(base_size = 14)
dev.off()

##################
## Analysing differences between causation and correlation
######################
meta.cor = read.table('../data_files/meta_cor.txt', sep = '\t', header = T)

pval_difference  =  2 * pnorm( -abs( (meta.cor$TE.fixed - meta.cor$Cor) / sqrt(meta.cor$seTE.fixed^2 + meta.cor$moyenne^2) ))
p = p.adjust(pval_difference, method = 'bonferroni')
meta.cor$t.test = pval_difference
meta.cor$t.test.adjust = p
write.table(meta.cor, '../data_files/meta_cor.txt', sep = '\t', col.names = T, row.names = F)

######################
## END
######################
print("===================")