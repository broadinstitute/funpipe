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
    version='0.1',
    description='A pipeline for analyzing fungal genomic data',
    long_description=readme(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Topic :: Bioinformatics',
    ],
    keywords=['fungal genomics processing pipeline'],
    url='https://github.com/broadinstitute/funpipe',
    project_urls={
        'Bug tracker': 'https://github.com/broadinstitute/funpipe/issues',
        'Documentation': 'https://github.com/broadinstitute/funpipe/README.rst'
    },
    author='Xiao Li',
    author_email='xiaoli.cbs@gmail.com',
    license='MIT',
    packages=['funpipe', 'funpipe.scripts'],
    package_dir={
        'funpipe': 'src',
        'funpipe.scripts': 'scripts'
    },
    install_requires=[
        'argparse', 'crimson', 'pandas', 'matplotlib'
    ],
    test_suite='setup.test_suite',
    scripts=glob('scripts/*'),
    include_package_data=True,
    zip_safe=False
)
