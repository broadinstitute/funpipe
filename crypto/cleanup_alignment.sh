##!/bin/sh
# clean up alignment pilot

wkdir='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/batch1_AD_coverage'
cd $wkdir

samtools view -bh -f 2 Y290-90.indels_realigned.sorted.bam > Y290-90.indes_realigned_proper_aligned_sorted.bam
samtools flagstat Y290-90.indes_realigned_proper_aligned_sorted.bam

samtools view Y290-90.indes_realigned_proper_aligned_sorted.bam | awk '$3 ~ /_A/ {print $0}' | less
nohup samtools view -h Y290-90.indes_realigned_proper_aligned_sorted.bam | awk '$3 ~ /_A/ {print $0}' > Y290-90.indes_realigned_proper_aligned_sorted_A.sam &
