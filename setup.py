#!/usr/bin/env python

# setup script for mmaeve

from distutils.core import setup
from codecs import open
from os import path, walk

setup(
    name="mmaeve",
    version="0.5",
    description="A module for creating biomolecular structures in " + \
                "Python",
    url="https://github.com/AppleIntusion/MMAEVe",
    author="Gubbin Eel",
    license="GPLv2",
    requires=["numpy", "scipy", "pandas"],
    packages=["mmaeve"]
)
