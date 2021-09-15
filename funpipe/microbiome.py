import os
import sys
sys.path.append('.')
from utils import run

def diamond_blastx(fa, output):
    ''' blastx with diamond
    
    Parameters
    ----------
    fa: string
        The path to input fasta file.
    output: string
        The path to output file.
    
    Returns
    -------
    string
        The path to output file.
        
    '''
    if not os.path.exists( fa ):
        raise Exception('Sorry, input fasta file does not exist')
    else:
        cmd = ' '.join(['diamond blastx -d nr -q', fa, '-o', output])
        run(cmd)
        
    return output
