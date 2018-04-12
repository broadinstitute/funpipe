# Cryptococcus neoformans serotype D project
Commands to reproduce analysis of Cryptococcus neoformans serotype D project.

## Introduction
Fungal genomic analysis pipeline, including pilon based reference genome refinement, snpeff variant annotation.

## Dependencies
* Analyses were set to run with `UGER` cluster. `dotkit` is needed to set proper tools in the running environment.  
* Python virtual environment is needed to run `widdler`.
* Set `fungal-pipeline/src` in your `PATH`.

### Reference genomes
Download reference genome for JEC21 from NCBI, perform md5sum check and patch contig names
```
get_annot.py              # in devel
python3 patch_chrs.py     # patch : need to add command line args
```

### Evaluate pilon improvement of JEC21 using illumina data
Reproduce the analysis, submit below jobs consecutively.
Working directory and reference files are:
```
dir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/  # working directory
qsub run_pilon.sh
```

To run the pilon pipeline manually, do:
```shell
python3 run_pilon.py \
  --outdir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21 \
  --fa /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/JEC21.fasta \
  --bam /seq/picard_aggregation/G138688/AFA_1003_15/v1/AFA_1003_15.bam \
  --prefix AFA_1003_15 \
  --threads 1 \
  --ram 16
```

To config snpEff database, run the following script.
```
qsub create_snpeff_db.sh
```
### Run snpEff

`snpEff v4.1g` was used for this task. snpEff is known to throw out pseudogenes and some transcripts on mitochondrion, and this has not been updated in latest version (v4.3t).
```
sh snpeff_db.sh
```

### GATK variant calling
Python virtual environment will be needed to use BTL's widdler. To setup the environment:
```
# create virtual environment
use Python-2.7
mkdir /cil/shed/sandboxes/user/ENV
cd /cil/shed/sandboxes/user
virtualenv ENV                          

# install dependent packages
source ENV/bin/activate
pip install ratelimit
pip install python-dateutil
pip install pytz
pip install sqlalchemy
```
Use tool `widdler` to automatic setup environment and analysis.

Launch GATK jobs:
```
sh run_gatk.sh gatk_single_sample.json # test run on single sample
sh run_gatk.sh gatk_pilot_samples.json # submit job
```
Monitor jobs:
```
widdler monitor         # monitor job
widdler query           # query all finished jobs
widdler abort <task id> # abort task
task_dir <task id>      # print task directory
```
Default task were sent to `gscid-cromwell`, if not working, could use alternative, such as `ale1`.
### Ploidy analysis
```
qsub run_ploidy.sh
```