FUNPIPE: a python library for building best practice fungal genomic analysis pipeline
-----
FUNPIPE is a python library, as well as several packaged tools to facilitate easy development of fungal genomic analysis pipelines. using state-of-the-art computational biology tools

## Requirements
* Python >= 3.4
* Pandas
* [crimson](https://github.com/bow/crimson): a library for parsing outputs of bioinformatics tools
* A list of bioinformatics tools need to be properly installed and export to `PATH`.

### Installation
To build functional pipelines with this library, virtual machine or virtual environment is recommended. Both `docker` and `conda` configs were provided to setup the proper environment:

```
# use docker

# conda

```
Note:
* `diamond=0.9.22` uses boost library, which depends on python2.7. This conflicts with funpipe's python version.
* `FreeTree` is not yet available in Bioconda. We are in the process of submitting it to the repo, and will keep it posted here.

Above two tools are not available if using the conda setup. Docker image will be created for `funpipe` to resolve this issue.

### Synposis
* [funpipe](./funpipe): a directory that contains python library
* [scripts](./scripts): a set of executables for high level analysis
* [tests](./tests): module tests
* `setup.py`: pip setup
* `conda_env.yml`: spec file for setting up conda environment

# Documentation
Below are major functionality of this pipeline, including
* Reference genome quality evaluation with `pilon`.
* Variant annotation using `snpEff`.
* Haploid variant calling `GATK`.
* Phylogenetic analysis.
* Ploidy analysis.

The following projects uses this library:
* [Cryptococcus neoformans serotype D project](https://github.com/broadinstitute/cneoformans_serod_analysis)
To use this library for your own research, please cite this manuscript.
### Evaluate reference genome quality with `pilon`: `run_pilon.py`

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

### Annotation genomic variants with `snpEff`: `run_snpeff.py`

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

### Run GATK-WDL
Before analysis, update `json` file to corresponding `gatk.wdl`, then launch analysis with:
```
widdler run gatk.wdl example.json
```

### Ploidy analysis
This analysis will only work with PERL version 5.14 - 5.18.
```
usage: run_ploidy.py [-h] -i BAM --faidx FAIDX [--prefix PREFIX] [-o OUT_DIR]

run ploidy analysis by analysing

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX       Prefix of output files
  -o OUT_DIR, --out_dir OUT_DIR
                        Output file

required arguments:
  -i BAM, --bam BAM     Input file
  --faidx FAIDX         fasta index
```

### Phylogenetic analysis
Phylogenetic analysis was performed firstly .
```

```

### Utility tools
There are several utility tools to facilitate the analysis, they locate under `utils` subdirectory.
`task_dir` returns the task directory of a specific on-prem submission.
```
task_dir <job ID>
```
`uges` will help login UGES queue.
`widdler` is a wrapper around `widdler.py` from Broad-BTL for easier utility. It automatically setup environment to run the job.
