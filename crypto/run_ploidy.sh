#!/usr/bin/env bash

#$ -l h_rt=02:00:00
#$ -N ploidy
#$ -l h_vmem=8g
#$ -pe smp 1
#$ -binding linear:1
#$ -R y
#$ -P gscid
#$ -t 1-162
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

# to do: input json file as input

# project directory
# prjdir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI
# prjdir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_ad
prjdir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI/batch1
prefix=batch1_86_D

# Input: bamlist
#bamlist=/cil/shed/sandboxes/xiaoli/fungal-pipeline/analysis/crypto/batch1_162_samples_bamlist.tsv
bamlist=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/test_ad/batch1_75_AD_realigned_bam_list.tsv

# faidx should match fasta files used for the alignment
faidx=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_H99_JEC21.fasta.fai
#faidx=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta.fai
#faidx=/seq/references/Cryptococcus_neoformans_JEC21/v0/Cryptococcus_neoformans_JEC21.fasta.fai
#faidx=/seq/references/Cryptococcus_neoformans_grubii_h99/v0/Cryptococcus_neoformans_grubii_h99.fasta.fai

# run analysis
bam=$(awk "NR==$SGE_TASK_ID" $bamlist | cut -f2)
echo "run_ploidy.py --bam $bam --out_dir $prjdir --faidx $faidx"
run_ploidy.py --bam $bam --out_dir $prjdir --faidx $faidx

cd $prjdir
# produce realigned bam list
paste <(ll *bam | awk '{print $8}' | cut -f1 -d.) <(ll *bam | awk '{print $8}') > "$prefix"_realigned_bam_list.tsv
# add scripts
