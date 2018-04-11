#!/usr/bin/env bash

source /broad/software/scripts/useuse
reuse .samtools-1.3
reuse Python-3.4
reuse .htslib-1.7

export PATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PATH
export PYTHONPATH=/cil/shed/sandboxes/xiaoli/fungal-pipeline/src:$PYTHONPATH

# patch JEC21
patch_ref_contigs.py \
  --dir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_JEC21 \
  --prefix GCF_000091045.1_ASM9104v1_genomic \
  --ctg_sufx _D

# patch H99
patch_ref_contigs.py \
  --dir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_H99 \
  --prefix GCF_000149245.1_CNA3_genomic \
  --ctg_sufx _A

# create AD fasta and gff files
cd /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly
cat NCBI_JEC21/GCF_000091045.1_ASM9104v1_genomic.patched.noMT.fasta NCBI_H99/GCF_000149245.1_CNA3_genomic.patched.noMT.fasta > NCBI_H99_JEC21.fasta
samtools faidx NCBI_H99_JEC21.fasta
