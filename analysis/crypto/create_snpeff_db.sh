#!/usr/bin/env bash

#$ -l h_rt=06:16:16
#$ -N snpeff_db
#$ -l h_vmem=16g
#$ -pe smp 1
#$ -binding linear:1
#$ -R y
#$ -P gscid
#$ -m ea
#$ -M xiaoli@broadinstitute.org
#$ -o /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db
#$ -wd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db
#$ -j y

source /broad/software/scripts/useuse
reuse Java-1.8

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH

snpeff_db.sh \
  . \
  /cil/shed/apps/external/snpEff/snpEff-4.1g \
  JEC21NCBI \
  /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.noMT.fasta \
  /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.noMT.gff \
  16
