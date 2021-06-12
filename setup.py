from setuptools import setup
from glob import glob
import unittest


def readme():
    with open('README.md') as f:
        return f.read()


def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite


setup(
    name='funpipe',
    version='0.1.0',
    description='A pipeline for analyzing fungal genomic data',
    long_description=readme(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English'
    ],
    keywords=['bioinformatics'],
    url='https://github.com/broadinstitute/funpipe',
    project_urls={
        'Bug tracker': 'https://github.com/broadinstitute/funpipe/issues',
        'Documentation': 'https://github.com/broadinstitute/funpipe/README.md'
    },
    author='Xiao Li',
    author_email='xiaoli@broadinstitute.org',
    license='MIT',
    packages=['funpipe', 'funpipe.scripts'],
    package_dir={
        'funpipe': 'funpipe',
        'funpipe.scripts': 'scripts'
    },
    install_requires=[
        'argparse>=1.1', 'crimson>=0.4.0', 'pandas>=0.23.4',
        'matplotlib>=3.0.2', 'seaborn>=0.9.0'
    ],
    test_suite='setup.test_suite',
    scripts=glob('scripts/*'),
    include_package_data=True,
    zip_safe=False
)
