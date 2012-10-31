#!/usr/bin/env python
# encoding: utf-8
from setuptools import setup
from sphinx.setup_command import BuildDoc


dependencies = """
numpy
pexpect
sphinx
"""

cmdclass = {'docs': BuildDoc}

setup(name='daopilot',
    version="0.0.1",
    author='Jonathan Sick',
    author_email='jonathansick@mac.com',
    description='DAOPHOT driver for python pipelines.',
    license='BSD',
    install_requires=dependencies.split(),
    cmdclass=cmdclass,
    packages=['daopilot'])
