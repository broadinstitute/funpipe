import os
import sys
import contextlib
from subprocess import check_call
import hashlib


    
def done(job_name):
    """ show that a job is done
    
    Parameters
    ----------
    job_name: string
        job name
        
    Returns
    -------
    int
        1
        
    """
    print(" - "+job_name+"is done.\n")
    return 1


@contextlib.contextmanager
def cd(dir):
    ''' change directory
    
    Parameters
    ----------
    dir: string
        new directory to change to
        
    '''
    original_path = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(original_path)


def run(cmd):
    ''' execute a specific command
    
    Parameters
    ----------
    cmd: string
        command to execute
        
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
    ''' remove a file
    
    Parameters
    ----------
    files: tuple of strings
        paths of files to remove
    
    Returns
    -------
    int
        1
        
    '''
    for file in files:
        run('rm '+file)
    return 1


def check_md5(file, checksum):
    """ perform md5 checksum to verify whether 2 files are the same.
    
    Parameters
    ----------
    file: string
        file to check
    checksum: string
        known checksum value in hexadecimal format
        
    Returns
    -------
    bool
        True if 2 files are the same, else False
        
    """
    with open(file,'rb') as file_to_check:
        data = file_to_check.read()    
        md5_returned = hashlib.md5(data).hexdigest()
        
    if md5_returned  == chechsum:
        return True
    else:
        return False
