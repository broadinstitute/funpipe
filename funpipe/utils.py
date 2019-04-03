import os
import sys
import contextlib
from subprocess import check_call


def done(job_name):
    print(" - "+job_name+"is done.\n")
    return 1


@contextlib.contextmanager
def cd(dir):
    ''' change directory
    :param dir: new directory to change to
    '''
    original_path = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(original_path)


def run(cmd):
    ''' execute a specific command
    :param cmd: command to execute
    '''
    sys.stderr.write(cmd+"\n")
    check_call(cmd, shell=True)
    print(" - Done: "+cmd+"\n")
    return 1


def rm(file):
    ''' remove a file
    :param file: path of file to remove
    '''
    for i in files:
        run('rm '+file)
    return 1


def check_md5(file, checksum):
    """ perform md5 check sum
    :param file: file to check
    :param checksum: known checksum value
    """
    md5 = check_call()
    if md5 == chechsum:
        return True
    else:
        return False
