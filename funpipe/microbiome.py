from .utils import run
"""
Diamond
=======
"""
def diamond_blastx(fa, output):
    ''' blastx with diamond
    
    Parameters
    ----------
    fa: string
        input fasta file.
    output: string
        output file
    
    Returns
    -------
    string
        output file
    '''
    cmd = ' '.join(['diamond blastx -d nr -q', fa, '-o', output])
    run(cmd)
    return output
