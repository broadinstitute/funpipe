library(argparse)

# To do:
#' Generalize is_D function in cov_plot
#' 


#' reverse coordinates
#' @description reverse coordinates to the reverse strand of 
#' @param vec input vector
#' @return a reverse vector of coordinates
#' Need to fix
#' 
reverse_coord <- function(vec) {
  n_elem <- length(vec)
  dist <- vec[n_elem] - vec[(n_elem-1):1]
  reverse_coord <- c(vec[1], vec[1] + dist)
  return(reverse_coord)
}

#' reorder position by defined chromosome orders. 
#' @param cov_df coverage data.frame
#' @param contig_order order of contigs
#' @return updated data.frame
#' 
reorder_chr <- function (cov_df, contig_order) {
  cov_df = cov_df[order(match(cov_df$chr, contig_order)),]
  cov_df$pos = 1:dim(cov_df)[1]
  return(cov_df)
}

#' get mid position of a position vector
mid_pos <- function(pos_vec) {
  return(mean(c(min(pos_vec), max(pos_vec))))
}

#' generate coverage_plot
#' @description 
#' @param cov coverage data.frame
#' @param sample sample name
#' @param col_vec vector of colors for each chromosome
#' @param contigs contigs need to plot out 
#' @param inverse_contig list of contigs need to reverse their physical postion
#' @param is_D: whether the subgenome is a D subgenome
#' 
cov_plot <- function (
  cov, sample, col_vec, contigs=NULL, inverse_contig=NULL, is_D=F) {
  cov1 <- cov[cov$chr %in% contigs, ]
  min_pos = min(cov1$pos)
  max_pos = max(cov1$pos)
  if(is_D) {
    plot(cov1$pos, cov1[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
         xlab=NA, ylab='cov', xlim=c(min_pos, max_pos))
  } else {
    plot(cov1$pos, cov1[[sample]], ylim=c(0,5), cex = 0, xaxt='n',
         xlab=NA, ylab='cov', main=sample, xlim=c(min_pos, max_pos))
  }
  # plot coverage for each contigs
  for (i in contigs) {
    if (is_D) {
      chr = i - 14
    } else {
      chr = i
    }
    if (!is.null(inverse_contigs) & (i %in% inverse_contigs)) {
      pos = rev(cov1$pos[cov1$chr==i])
      mtext(paste('chr', chr, '*', sep=''), side=1, 
            at=mid_pos(cov1$pos[cov1$chr==i]), cex=0.8)
    } else {
      pos = cov1$pos[cov1$chr==i]
      mtext(paste('chr', chr, sep=''), side=1, 
            at=mid_pos(cov1$pos[cov1$chr==i]), cex=0.8)
    }
    segments(pos, 0, pos, cov1[cov1$chr==i, sample], cols[i])
    grid(NA, 5, lwd=1.5)
  }
  return(1)
}

parser <- ArgumentParser(description='Generate coverage plot')
parser$add_argument('-d', '--density', help='input density tsv file')
parser$add_argument('-c', '--color', help='color file for each chromosomes')
parser$add_argument('-i', '--inverse_contigs', 
                    help='contigs that need to inverse positions to complment strand')
parser$add_argument('-p', '--prefix', help='output prefix')
args <- parser$parse_args()

# add option to process 
density <- 'batch1_AD_progeny.den'
density <- 'batch1_81_D.den'
density <- 'batch1_15AD_cov_fc.den'

density <- 'batch1_AD_progeny_cov_den.tsv'

prefix <- 'batch1_AD_progeny'
cov <- read.csv(density_file, sep='\t')
n_sample = dim(cov)[2] - 6

color <- read.csv('cneo_h99_jec21_chr_coord_5k.dat', sep='\t')
cols <- vector(length=28)
for (i in 1:dim(color)[1]) {
  cols[i] <- rgb(color$R[i], color$G[i], color$B[i])
}

# separate two sub-genomes
subg1 <- c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14)    # A 
subg2 <- c(15, 16, 25, 26, 18, 19, 20, 21, 23, 24, 17, 27, 28, 22)  # D 

cov = reorder_chr(cov, c(subg1, subg2))
cols <- vector(length=28)
color$chr = c(subg1, subg2)
for (i in 1:dim(color)[1]) {
  cols[color$chr[i]] <- rgb(color$R[i], color$G[i], color$B[i])
}

if (legacy) {
  start_column = 6
} else{
  start_column = 2
}

# inverse
inverse_contigs <- c(3, 4, 5, 7, 8, 10, 14)
png(paste0(prefix,'.cmp.png'), width=2500, height=n_sample*420, unit='px', 
    res=300)
par(mfrow=c(n_sample * 2, 1))
for (sample in names(cov)[(1+start_column):(n_sample+start_column)]) {
  print(sample)
  par(mar=c(1.5, 5, 1, 0.2))
  cov_plot(cov, subg1, sample, cols, inverse_contigs)
  cov_plot(cov, subg2, sample, cols, is_D=T)
}
dev.off()

# 
contigs = c(paste('chr', 1:14, '_A', sep=''), paste('chr', 1:14, '_D', sep=''))
pct_cov = read.csv('batch1_AD_progeny.pct_cov.tsv', sep='\t')
pct_cov = read.csv('batch1_81D_AD_cov_fc.pct_cov.tsv', sep='\t')
pct_cov = read.csv('batch1_15AD_cov.pct_cov.tsv', sep='\t')

prefix='batch1_15D'
pct_cov = pct_cov[order(match(pct_cov$contigs, contigs)),]
n_samples = dim(pct_cov)[2] - 1
png(paste(prefix,'cov.png', sep=''), width=3000, height=n_samples*500, 
    unit='px', res=300)
par(mfrow=c(n_samples, 2))
par(mar=c(2, 4, 2, 0.5))
for (sample in names(pct_cov)[2:n_samples]) {
  barplot(pct_cov[c(1:14), sample], names.arg=pct_cov$contigs[1:14],
          main=paste(sample, '_A', sep=''), ylab='% coverage', ylim=c(0,6))
  barplot(pct_cov[c(15:28), sample], names.arg=pct_cov$contigs[15:28],
          main = paste(sample, '_D', sep=''), ylab='% coverage', ylim=c(0,6))
}
dev.off()
