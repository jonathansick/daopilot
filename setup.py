#!/usr/bin/env python
# encoding: utf-8
from setuptools import setup
from sphinx.setup_command import BuildDoc


dependencies = """
numpy
pexpect
sphinx
"""

cmdclass = {'build_sphinx': BuildDoc}

setup(name='pydaophot',
    version="0.0.1",
    author='Jonathan Sick',
    author_email='jonathansick@mac.com',
    description='DAOPHOT interface for python',
    license='BSD',
    install_requires=dependencies.split(),
    cmdclass=cmdclass,
    packages=['pydaophot']
)
