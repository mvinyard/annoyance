import setuptools
import re
import os
import sys


setuptools.setup(
    name="annoyance",
    version="0.0.18",
    python_requires=">3.7.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="annoyance - single-cell AnnData wrapper of Spotify's Annoy library.",
    packages=setuptools.find_packages(),
    install_requires=[
        "annoy>=1.17.0",
        "anndata>=0.8",
        "pandas>=1.4.4",
        "pydk>=0.0.54",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
