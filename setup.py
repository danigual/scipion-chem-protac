"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-protac',  # Required
    version='0.1',  # Required
    description='Scipion plugin to model PROTAC-mediated ternary complexes.',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/danigual/scipion-chem-protac',  # Optional
    author='Daniel Gutiérrez',  # Optional
    author_email='danigual40@gmail.com',  # Optional
    keywords='scipion docking ppi protac scipion-3.0 cheminformatics',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'protac = protac'},
    package_data={  # Optional
       'protac': ['protocols.conf'],
    }
)
