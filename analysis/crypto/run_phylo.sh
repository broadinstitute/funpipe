#!/usr/bin/env bash

#$ -l h_rt=06:00:00
#$ -N phylo
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
reuse Python-2.7

PATH=/cil/shed/apps/internal/OBA/SNPs/Analysis:$PATH
PATH=/cil/shed/apps/external/phylogeny/FastTree-2.1.8:$PATH
export PATH=/cil/shed/apps/internal/OBA/SNPs:$PATH

prefix=crypto_serod_pilot.snp

cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI
cp /gsap/store/cromwell_executions/gatk/72f2438d-bda1-4ad8-9bcc-c941246141a6/call-FilterSnps/execution/filtered_SNPs.vcf ${prefix}.vcf

filterGatkGenotypes.py \
  --min_GQ 50 --min_percent_alt_in_AD 0.8 --min_total_DP 10 \
   ${prefix}.vcf > ${prefix}_GQ50_AD08_DP10.vcf

echo ${prefix}_GQ50_AD08_DP10.vcf > vcf_list.txt
vcfSnpsToFasta.py --max_amb_samples 10 vcf_list.txt > ${prefix}.fasta
FastTreeDP -nt ${prefix}.fasta > ${prefix}.nwk
