#!/usr/bin/env python

# setup script for mmaeve

from distutils.core import setup
from codecs import open
from os import path, walk

here = path.abspath(path.dirname(__file__))

setup(
    name="mmaeve",
    version="0.5",  # change this in nmrglue/__init__.py also
    description="A module for creating biomolecular structures in " + \
                "Python",
    url="https://github.com/AppleIntusion/MMAEVe",
    author="Gubbin Eel",
    license="New BSD License",
    requires=["numpy", "scipy", "pandas"],
    packages=["mmaeve"]
)
