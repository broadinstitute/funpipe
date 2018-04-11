#!/usr/bin/env bash

#$ -l h_rt=02:00:00
#$ -N ploidy
#$ -l h_vmem=8g
#$ -pe smp 1
#$ -binding linear:1
#$ -R y
#$ -P gscid
#$ -t 1-24
#$ -tc 0
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/ploidy
#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
reuse UGER
reuse .samtools-1.3
reuse Python-3.4
reuse .htslib-1.7
reuse .perl-5.18.1

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH
export PYTHONPATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PYTHONPATH

wkdir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454
prjdir="$wkdir"/analysis/JEC21_NCBI

bamlist="$prjdir"/pilot.realigned.bam_list.tsv
faidx="$wkdir"/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta.fai

bam=$(awk "NR==$SGE_TASK_ID" $bamlist)
run_ploidy.py --bam $bam --out_dir $prjdir/ploidy --faidx $faidx
