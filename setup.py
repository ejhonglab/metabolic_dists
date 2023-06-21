#!/usr/bin/env python3

from setuptools import setup


setup(
    name='metabolic_dists',
    version='0.1.0',
    py_modules=['reaction_distance'],
    # So far, shouldn't need install_requires=[...], as I think I'm only using stuff
    # installed in the metabolike conda environment by default.
)
