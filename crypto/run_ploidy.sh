#!/usr/bin/env bash

#$ -l h_rt=02:00:00
#$ -N ploidy
#$ -l h_vmem=8g
#$ -pe smp 1
#$ -binding linear:1
#$ -R y
#$ -P gscid
#$ -m ea
#$ -M xiaoli@broadinstitute.org
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI
#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
reuse UGER
reuse .samtools-1.3
reuse Python-3.4

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH
export PYTHONPATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PYTHONPATH

cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI

dir=ploidy
bamlist=pilot.aligned.bam_list.tsv
faidx=assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta.fai

while read bam; do
  run_ploidy.py --bam $bam --out_dir $dir --faidx $faidx
done < $bamlist
