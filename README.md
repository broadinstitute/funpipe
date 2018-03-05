# Evaluate pilon improvement of JEC21 using illumina data
Run pilon to evaluate the quality of JEC21

## Prerequisite
* Python3
* Access to UGER farm

## Usage
Reproduce this work from UGER farm
```
use UGER
qsub src/launch_to_uger.sh
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
