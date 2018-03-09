#!/usr/bin/env bash

#$ -l h_rt=06:16:16
#$ -N snpeff
#$ -l h_vmem=11g
#$ -pe smp 1
#$ -binding linear:4
#$ -R y
#$ -P gscid
#$ -m ea
#$ -M xiaoli@broadinstitute.org
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
#$ -j y
#$ -wd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21

source /broad/software/scripts/useuse
reuse UGER
reuse GCC-5.2
reuse Java-1.8

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH

cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
run_snpeff.py \
  -i AFA_1003_15.vcf.gz \
  -o AFA_1003_15.snpeff.vcf.gz \
  --ram 10 --genome_name JEC21 \
  -c /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db/snpEff.config \
  --jar /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db/snpEff.jar
