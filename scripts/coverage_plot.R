#!/usr/bin/env Rscript
library('argparse')

# To do:
#' Generalize is_D function in cov_plot
#'


#' @title reverse coordinates
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
#' @param cols vector of colors for each chromosome
#' @param contigs contigs need to plot out
#' @param inverse_contig list of contigs need to reverse their physical postion
#' @param is_D: whether the subgenome is a D subgenome
#'
cov_plot <- function (
  cov, sample, cols, contigs=NULL, inverse_contigs=NULL, is_D=F) {
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


filter_background <- function (cov){
  if (!all(c('sample_idx', 'cutoff') %in% names(cov))) {
    stop('sample_idx or cutoff not available in coverage list')
  }
  cov$df[,cov$sample_idx] = apply(
    cov$df[,cov$sample_idx], c(1,2), function(x) if(x< cov$cutoff) 0 else x)
  return(cov)
}


get_sample_idx <- function(cov) {
  cov$n_sample = dim(cov$df)[2] - cov$start_column
  cov$sample_idx = (1+cov$start_column):(cov$n_sample+cov$start_column)
  return(cov)
}


parse_cov<- function (density, legacy, cutoff) {
  cov = list()
  cov$cutoff <- cutoff
  if (legacy) {
    cov$start_column = 6
  } else{
    cov$start_column = 2
  }
  cov$df <- read.csv(density, sep='\t')
  cov <- get_sample_idx(cov)
  if (cov$cutoff > 0) {
    cov <- filter_background(cov)
  }
  return(cov)
}


coverage_plot<- function(cov, prefix, inverse_contigs, subg1, subg2, cols) {
  png(paste0(prefix,'.png'), width=2500, height=cov$n_sample*420, unit='px',
      res=300)
  par(mfrow=c(cov$n_sample * 2, 1))
  for (sample in names(cov$df)[cov$sample_idx]) {
    print(sample)
    par(mar=c(1.5, 5, 1, 0.2))
    cov_plot(cov$df, sample, cols, subg1, inverse_contigs)
    cov_plot(cov$df, sample, cols, subg2, is_D=T)
  }
  dev.off()
  return(1)
}


main <- function(density, prefix, color_csv, legacy, cutoff, inverse_contigs) {
  cov <-  parse_cov(density, legacy, cutoff)

  # separate two sub-genomes
  color <- read.csv(color_csv, sep='\t',
                    col.names=c('chr', 'start', 'end', 'R', 'G', 'B'), header=F
  )
  subg1 <- c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14)    # A
  subg2 <- c(15, 16, 25, 26, 18, 19, 20, 21, 23, 24, 17, 27, 28, 22)  # D

  cov$df = reorder_chr(cov$df, c(subg1, subg2))
  cols <- vector(length=28)
  color$chr = c(subg1, subg2)
  for (i in 1:dim(color)[1]) {
    cols[color$chr[i]] <- rgb(color$R[i], color$G[i], color$B[i])
  }
  coverage_plot(cov, prefix, inverse_contigs, subg1, subg2, cols)
}


# main
parser <- ArgumentParser(description='Generate coverage plot')
parser$add_argument('-i', '--cov_tsv', help='input density tsv file', 
                    required=T)
parser$add_argument('-c', '--color', help='color file for each chromosomes',
                    required=T)
parser$add_argument('-p', '--prefix', help='output prefix', required=T)
parser$add_argument('-l', '--legacy', help='legacy mode', action='store_true')
parser$add_argument('--cutoff', help='coverage cutoff to remove background', 
                    default=0)
parser$add_argument('--inverse_contigs', 
                    help='contigs that need to inverse positions to complementary strand',
                    default= c(3, 4, 5, 7, 8, 10, 14))
args <- parser$parse_args()
main(args$cov_tsv, args$prefix, args$color, args$legacy, args$cutoff,
     args$inverse_contigs)
