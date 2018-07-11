#!/usr/bin/env bash

#$ -N test
#$ -l h_rt=06:16:16
#$ -l h_vmem=8g
#$ -pe smp 1
#$ -binding linear:1
#$ -R y
#$ -q gscid
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
#$ -j y
#$ -wd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
#$ -l hostname=btl-c003

source /broad/software/scripts/useuse
reuse UGES
reuse Java-1.8
reuse Python-3.4

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH
echo "Done"
#cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21
#run_snpeff.py \
#  -i AFA_1003_15.vcf.gz \
#  -o AFA_1003_15.snpeff.vcf.gz \
#  --ram 80 --genome_name JEC21 \
#  -c /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db/snpEff.config  \
#  --snpeff_jar /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db/snpEff.jar
