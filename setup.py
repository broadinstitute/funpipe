from setuptools import setup


def readme():
    with open('README.rst') as f:
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
      url='https://github.com/broadinstitute/FunPipe',
      author='Xiao Li',
      author_email='xiaoli.cbs@gmail.com',
      license='MIT',
      packages=['funpipe'],
      install_requires=[
        'argparse', 'crimson', 'pandas', 'matplotlib'
      ],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
#      scripts=['scripts/*', 'utils/*'],
      include_package_data=True,
      zip_safe=False
)
