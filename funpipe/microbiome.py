from .utils import run
"""
Diamond
=======
"""
def diamond_blastx(fa, output):
    ''' blastx with diamond
    
    Parameters
    ----------
    arg1: string
        fa: input fasta file.
    arg2: string
        output: output file
    
    Returns
    -------
    string
        output file
    '''
    cmd = ' '.join(['diamond blastx -d nr -q', fa, '-o', output])
    run(cmd)
    return output
