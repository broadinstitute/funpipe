Biolego: a python library for efficient development of bioinformatic analysis pipelines
-----
## Introduction
**biolego** is a python library designed for efficient implementation of bioinformatics pipelines. It contains wrapper functions to popular tools, customized functions for common bioinformatics analysis, and command line tools developed using those functions. This package is developing initially to facilitate fungal genomic analysis, but most of the functions are generally applicable to other genomic analysis as well.

## Requirements
* Python >= 3.7
* Pandas >= 0.23.4
* Matplotlib >= 3.0.2
* [Crimson](https://github.com/bow/crimson) >= 0.4.0: a library for parsing outputs of bioinformatics tools
* A list of bioinformatics tools need to be properly installed and add to `PATH`. See here for the list and their versions `conda_env.yml`.

### Installation
We recommend to use conda to setup `biolego`'s dependencies. Make sure `conda` is available in your environment via `which conda`. If `conda` is not available in your system, install Python3.7 version of it [here](https://conda.io/miniconda.html). For Broad Internal users, you can use dotkit `use Anaconda3`.

Setup biolego working environment.
```sh
# work with the latest version of biolego
git clone git@github.com:broadinstitute/biolego.git

# setup environment
cd biolego
conda env create -f conda_env.yml  # this will take about 10 min
conda list  # verify new environment was installed correctly

# install biolego in the virtual environment
source activate biolego
pip install .

# deactivate the environment when done
source deactivate

# completely remove the virtual environment
conda remove -name biolego --all
```
Note:
* `diamond=0.9.22` uses boost library, which depends on python2.7. This conflicts with biolego's python version. Docker image will be created for `biolego` to resolve this issue.

### Synposis
* [biolego](./biolego): a directory that contains python library
* [scripts](./scripts): a set of executables for high level analysis
* [tests](./tests): module tests
* `setup.py`: pip setup
* `conda_env.yml`: spec file for setting up conda environment

### Documentation
Below are major functionality of this pipeline, including
* Reference genome quality evaluation with `pilon`.
* Variant annotation using `snpEff`.
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

The following projects uses this library:
* [*Cryptococcus neoformans* serotype D project](https://github.com/broadinstitute/fungal-research-projects/blob/master/docs/crypto_d.md)
* [*Candida auris* global project](https://github.com/broadinstitute/fungal-research-projects/blob/master/docs/cauris_global.md)
