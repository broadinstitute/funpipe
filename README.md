FunPipe: a python library for building best practice fungal genomic analysis pipeline
-----
`FunPipe` is a python library designed for efficient implementation of bioinformatic tools and pipelines for fungal genomic analysis. It contains wrapper functions to popular tools, customized functions for specific analyses tasks, and command line tools developed using those functions. This package is developing to facilitate fungal genomics, but many of the functions are generally applicable to other genomic analysis as well.

## Requirements
* Python >= 3.7
* Pandas >= 0.23.4
* Matplotlib >= 3.0.2
* [Crimson](https://github.com/bow/crimson) >= 0.4.0: a library for parsing outputs of bioinformatics tools
* [Conda](https://conda.io/miniconda.html)
* [Bioinformatic tool collections](./conda_env.yml)

The above list of bioinformatic tools need to be properly installed and add to `PATH`. See `conda_env.yml` for the list and their versions.

### Installation
**Install with PIP**
Part of PIP install process will use [`conda`](https://conda.io) for the bioinformatic tool collections.
Make sure `conda` is available in your environment via `which conda`. If `conda` is not available in your system, install Python3.7 version of it [here](https://conda.io/miniconda.html).

```sh
pip install funpipe

# activate conda environment
conda activate funpipe

# deactivate the environment when done
conda deactivate
```
Note:
* `diamond=0.9.22` uses boost library, which depends on `python 2.7`. This conflicts with `funpipe`'s python version. Use dimond via docker.

**Setup via Conda**
To use the latest version of funpipe, you can set it up via `conda`.

```sh
# clone this repo
git clone git@github.com:broadinstitute/funpipe.git

# setup environment
cd funpipe
conda env create -f conda_env.yml  # this will take about 10 min
conda list  # verify new environment was installed correctly

# install funpipe in the virtual environment
conda activate funpipe
pip install .   # to do: need to avoid conda installation again

# deactivate the environment when done
conda deactivate

# completely remove the virtual environment
conda remove -name funpipe --all
```
Note:
* `diamond=0.9.22` uses boost library, which depends on `python 2.7`. This conflicts with funpipe's python version. To use diamond, use it via docker.

**Setup via Docker**
There's a bit more overhead using Docker, but it came along with the benefits of consistent analysis environment (including the operation systems). It's extremely useful when using `funpipe` on the cloud.

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

### Synposis
* [funpipe](./funpipe): a directory that contains python library
* [scripts](./scripts): a set of executables for high level analysis
* [tests](./tests): module tests
* `setup.py`: pip setup script
* `conda_env.yml`: spec file for setting up conda environment
* `Dockerfile`: docker images

### Documentation
Major analysis pipelines:
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
    - Structural variant analysis with `breakdancer`
- Phylogenetic analysis
  - Dating analysis with `BEAST`.
  - Phylogenetic tree with `FastTree`, `RAxML` and `IQTREE`.
- GWAS analysis with `GEMMA`.

Here are scripts to run each of the above pipelines, use `<toolname> -h` to see manual.
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
You can also use out APIs to build your customized analysis scripts or pipelines. Checkout documents at: https://funpipe.readthedocs.io
