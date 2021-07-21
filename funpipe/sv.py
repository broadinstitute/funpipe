from .utils import run
"""
Breakdancer
===========
"""

def breakdancer(bam_file, prefix):
    ''' Detect structural variation using breakdancer
    
    Parameters
    ----------
    arg1: string
        bam_file: input bam file
    arg2: string
        prefix: output prefix
        
    Returns
    -------
    string
        contig file
    '''
    # create config files
    cfg_file = prefix+'.cfg'
    run('bam2cfg.pl -g -h'+bam_file+'> '+prefix+'.cfg')
    # Detect chromosomal structural variants using breakdancer-max
    run('brakdancer-max -q 40 -r 20 -y 90 '+cfg_file)
    return
