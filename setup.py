#!/usr/bin/env python

from distutils.core import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='ionolib',
      version='0.1',
      description='Slice, dice, and visualize numerical ionospheres',
      author='Nathaniel A. Frissell',
      author_email='nathaniel.frissell@scranton.edu',
      url='https://hamsci.org',
      packages=['ionolib'],
      install_requires=requirements
     )
