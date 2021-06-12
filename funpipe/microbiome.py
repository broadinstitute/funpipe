from .utils import run


def diamond_blastx(fa, output):
    ''' blastx with diamond
    :param fa: input fasta file.
    :param output: output file
    :return output
    '''
    cmd = ' '.join(['diamond blastx -d nr -q', fa, '-o', output])
    run(cmd)
    return output
