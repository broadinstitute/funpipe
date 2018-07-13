##!/bin/sh
# clean up alignment pilot

wkdir='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/batch1_AD_coverage'
cd $wkdir

samtools view -bh -f 2 Y290-90.indels_realigned.sorted.bam > Y290-90.indes_realigned_proper_aligned_sorted.bam


#
samtools view -bh Y290-90.indels_realigned.sorted.bam chr1_A > Y290-90.indels_realigned.sorted.chr1_A.bam
samtools index Y290-90.indels_realigned.sorted.chr1_A.bam
samtools bedcov chr1_A.bed Y290-90.indels_realigned.sorted.chr1_A.bam > Y290-90.indels_realigned.sorted.chr1_A.txt

#
samtools view -bh Y290-90.indes_realigned_proper_aligned_sorted.bam chr1_A > Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam
samtools index Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam
samtools bedcov chr1_A.bed Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam > Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.txt

samtools view -bh -f 2 -F 1024 -F 4 -F 8 -F 512 -F 2048 -F 256 -q 30 Y290-90.indels_realigned.sorted.bam chr1_A > Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam
samtools index Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam
samtools bedcov chr1_A.bed Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam > Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.txt
# explore high coverage regions
samtools view -bh Y290-90.indes_realigned_proper_aligned_sorted_chr1_A.bam chr1_A:10000-15000 > Y290-90.indes_realigned_proper_aligned_sorted_chr1_A_10000.bam


samtools flagstat Y290-90.indes_realigned_proper_aligned_sorted.bam

samtools view Y290-90.indes_realigned_proper_aligned_sorted.bam | awk '$3 ~ /_A/ {print $0}' | less
nohup samtools view -h Y290-90.indes_realigned_proper_aligned_sorted.bam | awk '$3 ~ /_A/ {print $0}' > Y290-90.indes_realigned_proper_aligned_sorted_A.sam &


# R plot:

setwd('~/Desktop')
cov = read.csv('290-90.indes_realigned_proper_aligned_sorted_chr1_A.txt', sep='\t', header=F)
cov$ave_cov = cov$V4/5000
plot(cov$V2, cov$ave_cov, xlab = 'pos', ylab = 'Average Coverage')
segments(cov$V2, 0, cov$V2, cov$ave_cov)
