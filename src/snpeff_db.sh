#!/bin/bash
set -e

dir=$1             # destination for snpdb
snpeff=$2          # path for snpeff jar
genome_name=$3     # name of genome
ref_fa=$4          # path to fasta file
gff3=$5            # path to gff3 file
ram=$6             # RAM usage

cd $dir
rm -rf *
# snpeff='/cil/shed/apps/external/snpEff/snpEff-4.3t/snpEff'
# snpeff='/cil/shed/apps/external/snpEff/snpEff-4.1g'
# snpeff='/seq/annotation/bio_tools/snpEff_2_0_5'

# genome_name='JEC21' # name of reference
# ref_fa='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.fna'      # reference genome
# gff3='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.gff'

# Please update the config file by:
sed 's/genomes : \\/genomes : JEC21, \\/g' $snpeff"/snpEff.config" > snpEff.config
echo $genome_name'.genome: Cryptococcus_neoformans_serotype_D' >> snpEff.config

# need to manually update config file
# genomes : ailmel1.61, ERTm2, your_genome_name

# Need to create a symlink to the jar, otherwise, it will use config file from
# where it locates and re-write config file in the current directory.
if [ ! -e snpEff.jar ]; then
  ln -s $snpeff"/snpEff.jar" .
fi

# Construct snpEff database.
if [ ! -d data ]; then
  mkdir data
  mkdir data/genomes
fi

mkdir data/$genome_name

# create symlinks
cd data
ln -s $ref_fa genomes/$genome_name.fa
ln -s $gff3 $genome_name/genes.gff
cd ..

java -Xmx${ram}g -jar snpEff.jar build -gff3 -v $genome_name \
    -config snpEff.config > ${genome_name}_build.log 2>&1
echo "Done."
