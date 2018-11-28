# FunPipe: a python library for building best practice fungal genomic analysis pipeline

`FunPipe` is a python library designed for efficient implementation of bioinformatic tools and pipelines for fungal genomic analysis. It contains wrapper functions to popular tools, customized functions for specific analyses tasks, and command line tools developed using those functions. This package is developing to facilitate fungal genomics, but many of the functions are generally applicable to other genomic analysis as well.

## Synposis
* [funpipe](./funpipe): a directory that contains python library
* [scripts](./scripts): tools and established pipelines, doc [here](#DOC)
* [tests](./tests): unit tests
* [docs](./docs): API documentation
* `README.md`: this file
* `setup.py`: pip setup script
* `conda_env.yml`: spec file for setting up conda environment
* `Dockerfile`: docker images
* `requirements.txt`: sphinx requirement file (not requirement for this package)
* `LICENSE`: MIT license

## Installation
### **<a name='CONDA'>Install with Conda</a>**
It is recommended to install funpipe via `conda`, as it automatically setup all required bioinformatic tools. This is very useful on servers or
clusters without root privilage. Make sure `conda` is available in your environment via `which conda`. If `conda` is not available in your system, install Python3.7 version of it [here](https://conda.io/miniconda.html).

HTTP errors sometimes occur when creating the conda environment, simply rerun the `conda env create -f conda_env.yml` to continue creating the environment.

```sh
# clone this repo
git clone git@github.com:broadinstitute/funpipe.git

# setup conda environment
cd funpipe

conda env create -f conda_env.yml # this will take a few mins
conda list  # verify new environment was installed correctly

# activate funpipe environment
conda activate funpipe

# the latest stable version of funpipe is available in this environment
# to use the latest funpipe version, do
pip install .

# deactivate the environment when done
conda deactivate

# completely remove the virtual environment
conda remove -name funpipe --all

# to complete remove the environment
conda remove -n funpipe --all
```

Note:
* `diamond=0.9.22` uses boost library, which depends on `python 2.7`. This conflicts with funpipe's python version. To use diamond, use it via [docker](#DOCKER).

### **<a name='DOCKER'>Install via Docker</a>**
There's a bit more overhead using Docker, but it came along with the benefits of consistent  environment (i.e.: including the operation systems). It's very useful when using `funpipe` on the cloud.

To use docker:
```
# Download docker
docker pull broadinstitute/funpipe:latest

# Run analysis interactively
docker run --rm -v $path_to_data/data -t broadinstitute/funpipe \
    /bin/bash -c "/scripts/vcf_qc_metr.py \
        -p prefix --jar /bin/GenomeAnalysisTK.jar \
        --fa /data/reference.fa
    "
```

You can use `Dockerfile` to compile the docker from scratch:
```sh
cd funpipe
docker build funpipe .
```

### **Install with PIP**
This approach is for advanced users who don't like conda and want to integrate funpipe into their current working environment. Before starting pip installation, make sure the following list of bioinformatic tools (or a subset of tools of interest) are properly installed and add to your `PATH`. Path to Java tools (JARs) need to be specified when evocaking specific functions.

**Requirements**
* Python >= 3.7
* Bioinformatic tool collections: can be automatically installed via conda [here](#CONDA)
    * Basic functions:
        - samtools>=1.9
        - bwa>=0.7.8
        - gatk>=3.8
        - picard>=2.18.17
    * Phylogenetics:
        - raxml>=8.2.12
        - readseq>=2.1.30
    * CNV:
        - breakdancer>=1.4.5
        - cnvnator>=0.3
    * Microbiome:
        - pilon>=1.23
        - diamond>=0.9.22

To install with pip:
```sh
# install latest stable release
pip install funpipe

# install a specific version
pip install funpipe==0.1.0
```

To install the latest version: funpipe
```sh
git clone git@github.com:broadinstitute/funpipe.git
cd funpipe
pip install .
```

## <a name='DOC'>Documentation</a>
Major analysis pipelines/tools:
- Quality control modules
    - Reference genome quality evaluation with `Pilon`.
    - FASTQ quality control with `fastqc`.
    - BAM quality control using `Picard`.
    - VCF quality control using `GATK VariantEval`.
- Variant Annotation with `snpEff`.
- Genomic Variation
    - Coverage analysis
    - Mating type analysis
    - Copy number variation with `CNVnator`
- Phylogenetic analysis
  - Dating analysis with `BEAST`.
  - Phylogenetic tree with `FastTree`, `RAxML` and `IQTREE`.
- GWAS analysis with `GEMMA`.

Here are scripts to run each of the above pipelines, use `<toolname> -h` to see the manuals.
```sh
##### Quality control #####
run_pilon.py          # Evaluate reference genome quality with pilon
fastqc.py             # Fastq quality control
bam_qc_metr.py        # Quality control of BAMs
vcf_qc_metr.py        # Quality control of VCFs

##### Variant Annotation #####
run_snpeff.py         # Annotation genomic variants with snpEff
phylo_analysis.py     # Phylogenetic analysis

##### Genomic Variations #####
coverage_analysis.py  # Hybrid coverage and ploidy analysis

```
You can also use out APIs to build your customized analysis scripts or pipelines. The docs will be available here: https://funpipe.readthedocs.io
