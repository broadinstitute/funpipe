from setuptools import setup
from glob import glob

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='funpipe',
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
        "Bug tracker": "https://github.com/broadinstitute/funpipe/issues",
        "Documentation": "https://github.com/broadinstitute/funpipe/README.rst"
      },
      author='Xiao Li',
      author_email='xiaoli.cbs@gmail.com',
      license='MIT',
      packages=['funpipe'],
      package_dir = {"funpipe": "src"},
      install_requires=[
        'argparse', 'crimson', 'pandas', 'matplotlib'
      ],
      test_suite='tests',
      scripts=glob('scripts/*'),
      include_package_data=True,
      zip_safe=False
)
