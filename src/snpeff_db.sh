#!/bin/bash
set -e
dir='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/annotation/snpeff_db'
cd $dir
rm -rf *
# snpeff='/cil/shed/apps/external/snpEff/snpEff-4.3t/snpEff'
# snpeff='/cil/shed/apps/external/snpEff/snpEff-4.1g'
snpeff='/seq/annotation/bio_tools/snpEff_2_0_5'

genome_name='JEC21' # name of reference
ref_fa='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.fna'      # reference genome
gff3='/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/NCBI_05Mar18/GCF_000149245.1_CNA3_genomic.patched.gff'

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
ln -s $ref_fa genomes/JEC21.fa
ln -s $gff3 $genome_name/genes.gff
cd ..

java -Xmx8g -jar snpEff.jar build -gff3 -v JEC21 -config snpEff.config > JEC21_build_030718.log 2>&1
echo "Done."
