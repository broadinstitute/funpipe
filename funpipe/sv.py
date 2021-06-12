from .utils import run


def breakdancer(bam_file, prefix):
    ''' Detect structural variation using breakdancer
    :param bam_file: input bam file
    :param prefix: output prefix
    :return
    '''
    # create config files
    cfg_file = prefix+'.cfg'
    run('bam2cfg.pl -g -h'+bam_file+'> '+prefix+'.cfg')
    # Detect chromosomal structural variants using breakdancer-max
    run('brakdancer-max -q 40 -r 20 -y 90 '+cfg_file)
    return
