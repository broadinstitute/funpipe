#!/usr/bin/env Rscript
library('argparse')
library('hash')
library('covR')

# To do:
.is.cov <- function(obj) {
  if (class(obj) != 'cov' ){
    stop("Input should be an coverage object")
  }
}

#' Subgenome matching between contigs in A and D
remove_subg <- function (contig_vec) {
  contig_vec <- gsub('_[AD]', '', contig_vec)
  contig_vec <- gsub('chr', '', contig_vec)
  return(as.factor(as.numeric(contig_vec)))
}

subg_map <- function() {
  subg_a <- paste('chr',c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), '_A',
                  sep='')
  subg_d <- paste('chr',c(1, 2, 11, 12, 4, 5, 6, 7, 9, 10, 3, 13, 14, 8), '_D',
                  sep='')
  return(hash::hash(subg_d, subg_a))
}


reorder_contigs <- function(vec, hs){
  for (key in hash::keys(hs)){
    if (key %in% vec) {
      idx <- which(vec == key)
      vec[idx] <- hs[[key]]
    }
  }
  return(vec)
}


#' @title reverse coordinates
#' @description reverse coordinates to the reverse strand, for ease of
#' comparing two subgenomes with different strand order
#' @param vec input coordinate vector
#' @return a reverse vector of coordinates
#' @example
#' vec <- c(1, 3, 6, 10)
#' rvec <- reverse_coord(vec)
#' rvec
#' [1] 1 5 8 10
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
#' @example
reorder_chr <- function (cov_df, contig_order) {
  cov_df = cov_df[order(match(cov_df$chr, contig_order)), ]
  cov_df$pos = 1:dim(cov_df)[1]
  return(cov_df)
}

#' get mid position of a position vector
mid_pos <- function(pos_vec) {
  return(mean(c(min(pos_vec), max(pos_vec))))
}

#' generate coverage barplot
#' @description
#' @param cov_df coverage data.frame
#' @param sample sample name
#' @param cols vector of colors for each chromosome
#' @param contigs contigs need to plot out
#' @param inverse_contig list of contigs need to reverse their physical postion
#' @param is_D: whether the subgenome is a D subgenome
coverage_barplot <- function (
  cov_df, sample, cols, contigs=NULL, inverse_contigs=NULL, is_D=F) {
  cov1 <- cov_df[cov_df$chr %in% contigs, ]
  min_pos = min(cov1$pos)
  max_pos = max(cov1$pos)
  # generate plot frame
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


#' An example coverage dataframe
# cov_df_example <- function() {
#   cov_df <- data.frame(
#
#   )
#   return(cov_df)
# }


#' remove background noise
#' @param cov coverage list
cov.filter_background <- function(cov) {
  .is.cov(cov)
  if (!all(c('sample_idx', 'cutoff') %in% names(cov))) {
    stop('sample_idx or cutoff not available in coverage list')
  }
  cov$df[,cov$sample_idx] = apply(
    cov$df[,cov$sample_idx], c(1,2), function(x) if(x< cov$cutoff) 0 else x)
  return(cov)
}


#' @param cov coverage list
cov.get_sample_idx <- function(cov) {
  .is.cov(cov)
  cov$n_sample = dim(cov$df)[2] - cov$start_column
  cov$sample_idx = (1+cov$start_column):(cov$n_sample+cov$start_column)
  return(cov)
}


#' @return a coverage list
cov.new <- function (cov_tsv, cutoff) {
  .is.cov(cov)
  cov = list()
  cov$df <- read.csv(cov_tsv, sep='\t')
  cov$cutoff <- cutoff
  class(cov) <- 'cov'
  return(cov)
}

cov.get_contigs <- function(cov) {
  .is.cov(cov)
  return(unique(cov$df$chr))
}


#' Profile preprocessing of coverage data, including
#' a. clean up background signal
#' b. inverse positions of contigs on the reverse strand
#' c. reorder subgenome contigs for ease of comparsion
#' d. create
#' @param cov coverage list
#' @description Including reorder D subgenome, inverse to reverse strand
cov.preprocessing <- function (cov, cutoff, inverse_contigs) {
  .is.cov(cov)
  cov <- cov.get_sample_idx(cov)
  cov$cutoff <- cutoff
  if (cov$cutoff > 0) cov <- filter_background(cov)
  # whether from D subgenome
  cov$df$is_D <- cov$df$chr %in% SUBG2
  # create ploting contig
  cov$df$plt_order <-
  # inverse positions of
  if(!is.null(inverse_contigs)) {
    contigs <- cov.get_contigs(cov)
    for (i in contigs) {
      idx = which(cov$df$chr == i)
      if (!is.null(inverse_contigs) & (i %in% inverse_contigs)) {
        cov$df$pos[idx] = rev(cov$df$pos[idx])
      }
    }
  }
  # reorder chromosomes
  cov$df = reorder_chr(cov$df, c(SUBG1, SUBG2))
  return(cov)
}

#' Create barplot
cov.barplot<- function(cov, prefix, inverse_contigs, subg1, subg2, cols) {
  .is.cov(cov)
  options(bitmapType='cairo')
  png(paste0(prefix,'_coverage_scatterplot.png'), width=2500,
      height=cov$n_sample*420, unit='px', res=300)
  par(mfrow=c(cov$n_sample * 2, 1))
  for (sample in names(cov$df)[cov$sample_idx]) {
    par(mar=c(1.5, 5, 1, 0.2))
    coverage_barplot(cov$df, sample, cols, subg1, inverse_contigs)
    coverage_barplot(cov$df, sample, cols, subg2, is_D=T)
  }
  dev.off()
  return(1)
}


create_color_vector <- function(color_csv, subg1, subg2) {
  color <- read.csv(
    color_csv, sep='\t', col.names=c('chr', 'start', 'end', 'R', 'G', 'B'),
    header=F
  )
  cols <- vector(length=28)
  color$chr = c(subg1, subg2)
  for (i in 1:dim(color)[1]) {
    cols[color$chr[i]] <- rgb(color$R[i], color$G[i], color$B[i])
  }
  return(cols)
}


main <- function(cov_tsv, prefix, color_csv, cutoff, inverse_contigs) {
  cov <-  cov.new(cov_tsv)  # create coverage list object
  cov <- cov.preprocessing(cov, cutoff, inverse_contigs)  # preprocess cov list
  cols <- create_color_vector(color_csv, SUBG1, SUBG2)

  coverage_plot(cov, prefix, inverse_contigs, SUBG1, SUBG2, cols)
  coverage_overlay_plot(cov, prefix, inverse_contigs, SUBG1, SUBG2, cols)
}


# input parser
if (!interactive()) {
  parser <- ArgumentParser(description='Generate coverage plot')
  parser$add_argument('-i', '--cov_tsv', help='input coverage tsv file',
                      required=T)
  parser$add_argument('-c', '--color', help='color RGB code for each chromosomes',
                      required=T)
  parser$add_argument('-p', '--prefix', help='output prefix', required=T)
  parser$add_argument('-l', '--legacy', help='legacy mode', action='store_true')
  parser$add_argument('--cutoff', help='coverage cutoff to remove background',
                      default=0)
  parser$add_argument('--inverse_contigs',
                      help='contigs that need to inverse positions to complementary strand',
                      default= c(3, 4, 5, 7, 8, 10, 14))
  args <- parser$parse_args()
  main(args$cov_tsv, args$prefix, args$color, args$cutoff, args$inverse_contigs)
}
