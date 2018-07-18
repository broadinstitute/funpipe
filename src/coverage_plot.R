library(grid)
# add option to process 
density_file <- 'batch1_AD_progeny.den'
density_file <- 'batch1_81_D.den'
density_file <- 'batch1_15AD_cov_fc.den'
prefix <- 'batch1_AD_progeny'
cov <- read.csv(density_file, sep='\t')

color <- read.csv('cneo_h99_jec21_chr_coord_5k.dat', sep='\t')
cols <- vector(length=28)
for (i in 1:dim(color)[1]) {
  cols[i] <- rgb(color$R[i], color$G[i], color$B[i])
}
n_sample = dim(cov)[2] - 6


# regular plot
chrs <- unique(color$chr)
png(paste0(prefix,'.png'), width=4000, height=n_sample*200, unit = 'px', res=300)
par(mfrow=c(n_sample, 1))
for (sample in names(cov)[7:(n_sample+6)]) {
  print(sample)
  par(mar=c(0.2, 5, 0.2, 0.2))
  # generate plot
  plot(cov$pos, cov[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
       xlab=NA, ylab=sample)
  # plot coverage for each chromosome
  for (i in chrs) {
    pos = cov$pos[cov$chr==i]
    segments(pos, 0, pos, cov[cov$chr==i, sample], cols[i])
    grid(NA, 5, lwd=1.5)
  }
}
dev.off()

# separate two sub-genomes
subg1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)    # A 
subg2 <- c(15, 16, 25, 19, 20, 21, 22, 28, 23, 24, 17, 18, 26, 27)  # D 
png(paste0(prefix,'cmp.png'), width=2000, height=n_sample*250, unit='px', res=300)
par(mfrow=c(n_sample * 2, 1))
for (sample in names(cov)[7:(n_sample+6)]) {
  print(sample)
  par(mar=c(0.2, 5, 1, 0.2))
  # generate plot
  cov1 <- cov[cov$chr %in% subg1, ]
  plot(cov1$pos, cov1[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
       xlab=NA, ylab='Copy', main=sample, xlim=c(min(cov1$pos), max(cov1$pos)))
  # plot coverage for each chromosome
  for (i in subg1) {
    pos = cov$pos[cov$chr==i]
    segments(pos, 0, pos, cov[cov$chr==i, sample], cols[i])
    grid(NA, 5, lwd=1.5)
  }
  
  # 
  cov2 <- cov[cov$chr %in% subg2, ]
  plot(cov1$pos, cov1[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
       xlab=NA, ylab='Copy', , xlim=c(min(cov2$pos), max(cov2$pos)))
  for (i in subg2) {
    pos = cov$pos[cov$chr==i]
    segments(pos, 0, pos, cov[cov$chr==i, sample], cols[i])
    grid(NA, 5, lwd=1.5)
  }
}
dev.off()
