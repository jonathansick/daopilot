#!/usr/bin/env python
# encoding: utf-8
"""
Allstar class wrapper.

2012-05-05 - Created by Jonathan Sick
"""

import os
import sys
import pexpect


class Allstar(object):
    """Wrapper object for Peter Stetson's allstar program for doing psf
    photometry, given a psf model made by daophot.
    
    .. todo:: refactor to allow one to pass a daophot object and it uses the
       'last' objects to automatically make all the paths.
    
    .. note:: All the inputs and output paths should be in the same directory
       as the `inputImagePath`.
    """
    def __init__(self, inputImagePath, psfPath, apPhotPath, alsOutputPath,
            outputImagePath, shell="/bin/zsh", cmd="allstar"):
        super(Allstar, self).__init__()
        self.shell = shell
        self.cmd = cmd
        self.inputImagePath = inputImagePath
        self.psfPath = psfPath  # psf model path
        self.apPhotPath = apPhotPath  # aperture photometry input path
        self.alsOutputPath = alsOutputPath  # allstar output (photometry) path
        self.outputImagePath = outputImagePath  # star-subtracted image path
        
        # Delete old copies of the output files
        if os.path.exists(self.alsOutputPath):
            os.remove(self.alsOutputPath)
        if os.path.exists(self.outputImagePath):
            os.remove(self.outputImagePath)
        
        # Will be a pexpect instance running allstar
        self.allstar = None
    
    def run(self, timeout=30. * 60):
        """Runs an allstar session.

        :param timeout: time (seconds) to allow `allstar` to run before
           giving up.
        """
        # need to delete the .als file, otherwise allstar will ask
        # to overwrite it
        if os.path.exists(self.alsOutputPath):
            os.remove(self.alsOutputPath)
        
        self.allstar = pexpect.spawn('%s -c "cd %s;%s"' %
                (self.shell, self.cmd, os.path.dirname(self.inputImagePath)))
        self.allstar.logfile = sys.stdout  # DEBUG
        self.allstar.expect("OPT>")
        print self.allstar.before
        
        # TODO add options?
        
        self.allstar.sendline("")
        self.allstar.expect("Input image name:")
        self.allstar.sendline(os.path.basename(self.inputImagePath))
        # asks File with the PSF
        self.allstar.expect(":")
        self.allstar.sendline(os.path.basename(self.psfPath))
        # asks Input file .ap
        self.allstar.expect(":")
        self.allstar.sendline(os.path.basename(self.apPhotPath))
        # asks for .als output path
        self.allstar.expect(":")
        self.allstar.sendline(os.path.basename(self.alsOutputPath))
        # asks for path for output star-sub image
        self.allstar.expect(":")
        self.allstar.sendline(os.path.basename(self.outputImagePath))
        
        # wait up to 30 minutes for allstar to finish
        self.allstar.expect("Good bye.", timeout=timeout)
        # TODO, will this get rid of the allstar build-up?
        # self.allstar.sendcontrol('d')
        print "finished"
        self.allstar = None
