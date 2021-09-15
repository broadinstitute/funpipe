import os
import sys
import contextlib
from subprocess import check_call
import hashlib


    
def done(job_name):
    """ Show that a job is done.
    
    Parameters
    ----------
    job_name: string
        The job name.
        
    Returns
    -------
    int
        1
        
    """
    print(" - "+job_name+"is done.\n")
    return 1


@contextlib.contextmanager
def cd(dir):
    ''' Change working directory.
    
    Parameters
    ----------
    dir: string
        New directory to change to.
        
    '''
    original_path = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(original_path)


def run(cmd):
    ''' Execute a specific command and print out the command.
    
    Parameters
    ----------
    cmd: string
        The command to execute.
        
    Returns
    -------
    int
        1
        
    '''
    sys.stderr.write(cmd+"\n")
    check_call(cmd, shell=True)
    print(" - Done: "+cmd+"\n")
    return 1


def rm(*files):
    ''' Remove file(s).
    
    Parameters
    ----------
    files: tuple
        Path(s) of file(s) to remove.
    
    Returns
    -------
    int
        1
        
    '''
    for file in files:
        run('rm '+file)
    return 1


def check_md5(file, checksum):
    """ Perform md5 checksum to verify whether 2 files are the same.
    
    Parameters
    ----------
    file: string
        The file to check.
    checksum: string
        The known checksum value in hexadecimal format.
        
    Returns
    -------
    bool
        True if 2 files are the same, else False.
        
    """
    with open(file,'rb') as file_to_check:
        data = file_to_check.read()    
        md5_returned = hashlib.md5(data).hexdigest()
        
    if md5_returned  == chechsum:
        return True
    else:
        return False
