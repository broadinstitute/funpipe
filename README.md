# Fungal pipeline
Fungal genomic analysis pipeline
## Introduction
Fungal genomic analysis pipeline, including:
* Reference genome quality evaluation
* snpEff annotation
* Variant calling
* Phylogenetic tree

## Projects using this pipeline
* [Cryptococcus neoformans serotype D project](analysis/crypto/README.md)

## Dependencies
* Python-3.4
* Bioinformatics tools need to be in the environment `PATH`.

## Usage
Below are major functionality of this pipeline:
### Evaluate reference genome quality with pilon: `run_pilon.py`

```
usage: run_pilon.py [-h] --prefix PREFIX --bam BAM --fa FA [--outdir OUTDIR]
                    [--gff3 GFF3] [--ram RAM] [--threads THREADS]
                    [--picard_jar PICARD_JAR] [--snpeff_jar SNPEFF_JAR]
                    [--pilon_jar PILON_JAR] [--snpeff_db SNPEFF_DB]

optional arguments:
  -h, --help            show this help message and exit
  --gff3 GFF3           GFF3 annotation
  --ram RAM             RAM usage of input file
  --threads THREADS     Number of threads for pilon
  --picard_jar PICARD_JAR
                        Picard jar
  --snpeff_jar SNPEFF_JAR
                        Jar to snpeff
  --pilon_jar PILON_JAR
                        Pilon Jar
  --snpeff_db SNPEFF_DB
                        snpEff database

Required arguments:
  --prefix PREFIX       output prefix
  --bam BAM             input bam
  --fa FA               reference genome fasta file to evaluate
  --outdir OUTDIR       output directory

```

### Annotation genomic variants with snpEff: `run_snpeff.py`

```
usage: run_snpeff.py [-h] -i INPUT_VCF -o OUTPUT_VCF --genome_name GENOME_NAME
                     -c CONFIG [-m RAM] [--snpeff_jar SNPEFF_JAR]

Run snpeff

optional arguments:
  -h, --help            show this help message and exit
  -m RAM, --ram RAM     RAM usage
  --snpeff_jar SNPEFF_JAR
                        jar file of snpeff

required arguments:
  -i INPUT_VCF, --input_vcf INPUT_VCF
                        input vcf
  -o OUTPUT_VCF, --output_vcf OUTPUT_VCF
                        Output vcf
  --genome_name GENOME_NAME
                        genome name in snpeff config file
  -c CONFIG, --config CONFIG
                        config file for snpeff
```

### WIDDLER wrapper
This wrapper automatically setup virtual environment for running WDL. 
```
widdler monitor         # monitor job
widdler query           # query all finished jobs
widdler abort <task id> # abort task
task_dir <task id>      # print task directory
```
