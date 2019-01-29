from . import utils


class microbiome(analysis):
    def __init__(self):
        analysis._init__()


def diamond_blastx(fa, output):
    """ blastx with diamond
    :param fa: input fasta file.
    :param output: output file
    :return output
    """
    cmd = ' '.join(['diamond blastx -d nr -q', fa, '-o', output])
    utils.run(cmd)
    return output
