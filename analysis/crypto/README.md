# Cryptococcus neoformans serotype D project
Documentation of analysis on Cryptococcus neoformans serotype D project.

## Introduction
Fungal genomic analysis pipeline, including pilon based reference genome refinement, snpeff variant annotation.

## Dependencies
Pipeline now can only be run on Broad's `UGER` cluster. Tools and environment are set using `dotkit`.  

### Evaluate pilon improvement of JEC21 using illumina data
Reproduce the analysis, submit below jobs consecutively (note, in the future will be a single shell script).
Working directory and reference files are:
```
dir=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/  # working directory
fa=/gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly
bams=
```

```
use UGER
qsub src/launch_to_uger.sh
qsub src/
qsub src/
```

To run the pilon pipeline manually, do:
```
python3 run_pilon.py \
  --outdir /gsap/garage-fungal/Crypto_neoformans_seroD_B454/analysis/JEC21 \
  --fa /gsap/garage-fungal/Crypto_neoformans_seroD_B454/assembly/JEC21.fasta \
  --bam /seq/picard_aggregation/G138688/AFA_1003_15/v1/AFA_1003_15.bam \
  --prefix AFA_1003_15 \
  --threads 1 \
  --ram 16
```

To config snpEff database, run the following script. Now that paths are hard-coded. SNPEff v2_0_5 was used for this task. snpEff is known to throw out pseudogenes and some transcripts, and this has not been updated in latest version (v4.3t) of it.
```
sh snpeff_db.sh
```
