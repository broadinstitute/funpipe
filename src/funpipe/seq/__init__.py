'''
Examples
--------
# use function
>>> from funpipe.seq.bam import sort_bam, index_bam
>>> soretd_sample = sort_bam(sample)

# use Bam Class
>>> from funpipe.seq.bam import Bam
>>> sample = Bam( 'sample.bam' )
Sort BAM file:
>>> sorted_sample = sample.sort_bam()
Clean BAM file:
>>> cleanup_sorted_sample = sorted_bam.clean_bam()
'''