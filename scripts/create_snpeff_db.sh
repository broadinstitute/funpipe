#!/bin/bash
set -o pipefail

crdir () {
    # create a directori if not exist
    if [ -d $1 ]; then
        echo "$1 already exist."
    else
        echo "creating $1"
        mkdir $1
    fi
}


usage () {
  cat <<HELP

  Create SNPEFF database from a FASTA and GFF file

  Usage:

  Use on UGER Queue via task array:
    snpeff_db.sh $dir $jar $genome_name $ref_fa $gff3 $ram
    dir: working directory
    jar: snpeff jar
    genome_name: name of the species
    ref_fa: path to the fasta file
    ram: int, ram allocated for the java VM.

  Note:
    create proper directory structure for snpeff to work
    can either refer to the absolute path of snpeff jar and specify config path,
    or create a symlink to the jar and run it from the local directory.

    ######## Structure ########
    # snpeff
    #  - snpEff.jar
    #  - snpEff.config
    #  - data
    #     - genomes
    #        - genome_name.fa
    #     - genome_name
    #        - genes.gtf
    ###########################

HELP
}


# main
dir=$1             # destination for snpdb
jar=$2             # path for snpeff jar
genome_name=$3     # name of genome
ref_fa=$4          # path to fasta file
gff3=$5            # path to gff3 file
ram=$6             # RAM usage

[[ -z $dir ]] && usage && exit 1

crdir ${dir}/snpeff
cd ${dir}/snpeff

# config
effdir=$(dirname $jar)
if [ -e $effdir/snpEff.config ]; then
    cp ${effdir}/snpEff.config .
else
    echo "snpeff config not available."
    exit 1
fi

# Construct snpEff database.
crdir data
crdir data/genomes
crdir data/$genome_name

ln -s $ref_fa data/genomes/${genome_name}.fa    # must be end with fa
ln -s $gff3 data/$genome_name/genes.gff

# Please update the config file by:

echo $genome_name".genome: $genome_name" >> snpEff.config

# need to manually update config file
# genomes : ailmel1.61, ERTm2, your_genome_name

# create database
java -Xmx${ram}g -jar $jar build -gff3 -v $genome_name \
    -config snpEff.config > ${genome_name}_build.log 2>&1

echo "Done."
