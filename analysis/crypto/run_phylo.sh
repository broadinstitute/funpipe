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
export $PATH

cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21_NCBI

prefix=crypto_serod_pilot
echo "${prefix}.vcf" > vcf_list.txt
vcfSnpsToFasta.py vcf_list.txt > ${prefix}.fasta
FastTreeDP -nt ${prefix}.fasta > ${prefix}.nwk
