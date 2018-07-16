library(grid)
density_file <- 'batch1_AD_progeny.den'
density_file <- 'batch1_81_D.den'
density_file <- 'batch1_15AD_cov_fc.den'
prefix <- 'batch1_15AD'
cov <- read.csv(density_file, sep='\t')

color <- read.csv('cneo_h99_jec21_chr_coord_5k.dat', sep='\t')
cols <- vector(length=28)
for (i in 1:dim(color)[1]) {
  cols[i] <- rgb(color$R[i], color$G[i], color$B[i])
}
chrs <- unique(color$chr)

# 
n_sample = dim(cov)[2] - 6
png(paste0(prefix,'.png'), width=4000, height=n_sample*200, unit = 'px', res=300)
par(mfrow=c(n_sample, 1))
for (sample in names(cov)[7:(n_sample+6)]) {
  print(sample)
  par(mar=c(0.2,5,0.2,0.2))
  plot(cov$pos, cov[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
       xlab=NA, ylab=sample)
  for (i in chrs) {
    pos = cov$pos[cov$chr==i]
    segments(pos, 0, pos, cov[cov$chr==i, sample], cols[i])
    grid(NA, 5, lwd=1.5)
  }
}
dev.off()

