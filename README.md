FUNPIPE: a python library for building best practice fungal genomic analysis pipeline
-----
## Introduction
**funpipe** is a python library designed for efficient implementation of bioinformatics pipelines. It contains wrapper functions to popular tools, customized functions for common bioinformatics analysis, and command line tools developed using those functions. This package is developing initially to facilitate fungal genomic analysis, but most of the functions are generally applicable to other genomic analysis as well.

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

# deactivate the environment when done
source deactivate

# completely remove the virtual environment
conda remove -name funpipe --all
```
Note:
* `diamond=0.9.22` uses boost library, which depends on python2.7. This conflicts with funpipe's python version. Docker image will be created for `funpipe` to resolve this issue.

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
* Coverage and ploidy analysis.

List of available tools, use `<toolname> -h` to see manual
```sh
run_pilon.py          # Evaluate reference genome quality with pilon
run_snpeff.py         # Annotation genomic variants with snpEff
fastqc.py             # Fastq quality control
phylo_analysis.py     # Phylogenetic analysis
coverage_analysis.py  # Hybrid coverage and ploidy analysis
bam_qc_metr.py        # Quality control of BAMs
vcf_qc_metr.py        # Quality control of VCFs
```
usage: run_pilon.py [-h] --prefix PREFIX --bam BAM --fa FA [--outdir OUTDIR]
                    [--gff3 GFF3] [--ram RAM] [--threads THREADS]
                    [--picard_jar PICARD_JAR] [--snpeff_jar SNPEFF_JAR]
                    [--pilon_jar PILON_JAR] [--snpeff_db SNPEFF_DB]

The following projects uses this library:
* [*Cryptococcus neoformans* serotype D project](https://github.com/broadinstitute/fungal-research-projects/blob/master/docs/crypto_d.md)
* [*Candida auris* global project](https://github.com/broadinstitute/fungal-research-projects/blob/master/docs/cauris_global.md)
