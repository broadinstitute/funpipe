#!/usr/bin/env bash

source /broad/software/scripts/useuse
use Python-2.7
source /cil/shed/sandboxes/xiaoli/ENV/bin/activate # activate virtual environment
export PATH=/cil/shed/sandboxes/xiaoli/widdler:$PATH

widdler.py run \
  /cil/shed/sandboxes/xiaoli/BTL-wdl/gatk/gatk.wdl gatk_pilot_samples.json \
  -S gscid-cromwell
