library(ggplot2)
library(grid)
library(gridExtra)
library(samplyzer)

qc = read.csv('~/Downloads/crypto_sero_d_BAM_D_to_D.tsv', sep='\t')
qc$index = 1:dim(qc)[1]

plts <- list()


plts[[1]] <- qplot(data=qc, index, MEAN_COVERAGE, color=Type)
plts[[2]] <- qplot(data=qc, index, PCT_PF_READS_ALIGNED, color=Type)
plts[[3]] <- qplot(data=qc, index, PCT_CHIMERAS, color=Type)
plts[[4]] <- qplot(data=qc, index, filtered_GQ, color=Type)
plts[[5]] <- qplot(data=qc, index, filtered_AD, color=Type)
plts[[6]] <- qplot(data=qc, index, filtered_DP, color=Type)

multiplotWithSharedLegend(plts, ncols=3, position = c('right'), show=T)

qc[qc$PCT_PF_READS_ALIGNED < 0.7,]
