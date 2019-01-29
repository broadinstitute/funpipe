from setuptools import setup
from glob import glob
import unittest
from os import path


def readme():
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
    return long_description


def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite


setup(
    name='biolego',
    version='0.0.1',
    description='a python library for efficient development of bioinformatic analysis pipelines',
    long_description=readme(),
    long_description_content_type='text/markdown',
    keywords = ['bioinformatics'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
    ],
    url='https://github.com/broadinstitute/biolego',
    project_urls={
        'Bug tracker': 'https://github.com/broadinstitute/biolego/issues',
        'Documentation': 'https://github.com/broadinstitute/biolego/README.md'
    },
    author='Xiao Li',
    author_email='xiaoli.cbs@gmail.com',
    license='MIT',
    packages=['biolego', 'biolego.scripts'],
    package_dir={
        'biolego': 'biolego',
        'biolego.scripts': 'scripts'
    },
    install_requires=[
        'argparse', 'crimson', 'pandas', 'matplotlib',
    ],
    test_suite='setup.test_suite',
    scripts=glob('scripts/*'),
    include_package_data=True,
    zip_safe=False
)
