#!/usr/bin/env bash

#$ -l h_rt=06:16:16
#$ -N JEC21_pilon
#$ -l h_vmem=16g
#$ -pe smp 1
#$ -binding linear:4
#$ -R y
#$ -P gscid
#$ -m ea
#$ -M xiaoli@broadinstitute.org
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
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
reuse GAEMER

export PATH=/gsap/garage-fungal/OBA/FungalProjects/B454_Crypto_PhaseIII/src:$PATH

run_pilon.py \
  --outdir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21 \
  --fa /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/JEC21.fasta \
  --bam /seq/picard_aggregation/G138688/AFA_1003_15/v1/AFA_1003_15.bam \
  --prefix AFA_1003_15 \
  --threads 1 \
  --ram 16
