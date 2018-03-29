#!/usr/bin/env bash

#$ -l h_rt=08:00:00
#$ -N filter_vcf
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
reuse Python-2.7

export PATH=/cil/shed/apps/internal/OBA/SNPs:$PATH

cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI
vcf=crypto_serod_pilot.vcf
cp /gsap/store/cromwell_executions/gatk/72f2438d-bda1-4ad8-9bcc-c941246141a6/call-CombineVariants/execution/filtered.combined.vcf $vcf

filterGatkGenotypes.py \
  --min_GQ 50 --min_percent_alt_in_AD 0.8 --min_total_DP 10 \
   $vcf > crypto_serod_pilot_GQ50_AD08_DP10.vcf

bgzip crypto_serod_pilot_GQ50_AD08_DP10.vcf
tabix crypto_serod_pilot_GQ50_AD08_DP10.vcf.gz
