# -*- coding: latin-1 -*-

import ez_setup
import sys
from setuptools import setup, find_packages

import elongwheat

"""
Notes:

- use setup.py develop when tracking in-development code
- when removing modules or data files from the project, run setup.py clean --all and delete any obsolete .pyc or .pyo.

"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

ez_setup.use_setuptools()

if sys.version_info < (2, 7):
    print('ERROR: Elong-Wheat requires at least Python 2.7 to run.')
    sys.exit(1)

if sys.version_info >= (3, 0):
    print('WARNING: Elong-Wheat has not been tested with Python 3.')

setup(
    name="Elong-Wheat",
    version=elongwheat.__version__,
    packages=find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'pandas'],

    # metadata for upload to PyPI
    author="M.Gauthier, C.Chambon, R.Barillot",
    author_email="camille.chambon@inra.fr, romain.barillot@inra.fr",
    description="Model of leaf elongation for wheat",
    long_description="A mechanistic model of leaf elongation for wheat that accounts for the CN status",
    license="CeCILL-C",
    keywords="functional-structural plant model, wheat, leaf growth, morphogenesis, trophic status, carbon, nitrogen, metabolism",
    url="https://sourcesup.renater.fr/projects/elong-wheat/",
    download_url="https://sourcesup.renater.fr/frs/download.php/file/5704/elongwheat.zip"
)
