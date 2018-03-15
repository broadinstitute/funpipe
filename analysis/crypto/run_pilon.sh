#!/usr/bin/env bash

#$ -l h_rt=06:16:16
#$ -N cserod_pilon
#$ -l h_vmem=16g
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
reuse GCC-5.2
reuse .samtools-1.3
reuse .bwa-0.7.12
reuse Java-1.8
reuse .python-3.6.3
reuse group=gtba
reuse GAEMR

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH

run_pilon.py \
  --outdir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI \
  --fa /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.noMT.fasta \
  --bam /seq/picard_aggregation/G138688/AFA_1003_15/v1/AFA_1003_15.bam \
  --prefix AFA_1003_15_NCBI \
  --threads 1 \
  --ram 16
