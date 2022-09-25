# funpipe

Calculate word counts in a text file!

## Installation

```bash
$ pip install funpipe
```

## Usage

`funpipe` can be used to count words in a text file and plot results
as follows:

### bam

```python
from funpipe.bam import bam

bam = bam( 'sample.bam' ) # path to BAM file
# Sort BAM file:
sorted_bam = bam.sort_bam().sorted_bam
# Clean BAM file:
cleanup_sorted_bam = sorted_bam.clean_bam().cleanup_bam
# Index BAM file:
cleanup_sorted_bam.index_bam()
# Compute the depth of alignment:
cleanup_sorted_bam.bam_depth( out_prefix = 'sample' )
# Detect structural variation:
cleanup_sorted_bam.breakdancer()
# Variant calling:
cleanup_sorted_bam.variant_calling( 'ref.fasta', 'sample')
```

## Contributing

Interested in contributing? Check out the contributing guidelines. 
Please note that this project is released with a Code of Conduct. 
By contributing to this project, you agree to abide by its terms.

## License

`funpipe` was created by x-lab. It is licensed under the terms
of the MIT license.

## Credits

`funpipe` was created with 
[`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and 
the `py-pkgs-cookiecutter` 
[template](https://github.com/py-pkgs/py-pkgs-cookiecutter).