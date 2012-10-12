#!/usr/bin/env python
# encoding: utf-8
"""
Class wrapper to daophot.

2012-05-05 - Created by Jonathan Sick
"""

import os
import sys

import pexpect


class Daophot(object):
    """Object-oriented interface to drive daophot.
    
    :param inputImagePath: is the path to the FITS image that will be measured.
        This is a real (filesystem) path. All paths should be supplied, and
        will be returned to the user as filesystem paths. The class internally
        converts these into shortened (symlinked) paths.
    """
    def __init__(self, inputImagePath):
        super(Daophot, self).__init__()
        self.inputImagePath = inputImagePath
        self.workDir = os.path.dirname(self.inputImagePath)
        # a pexpect process running daophot, None when shutdown
        self.daophot = None
        
        # Cache for paths; two levels of dictionaries. First level is keyed
        # to the types of files (represented by file extension strings). Second
        # level is keyed by the path names themselves
        self.pathCache = {'fits': {},
                'coo': {}, 'lst': {}, 'ap': {}, 'psf': {}, 'nei': {}}
        self.pathCache['fits']['input_image'] \
                = os.path.basename(self.inputImagePath)
        self.pathCache['fits']['last'] \
                = os.path.basename(self.inputImagePath)
    
    def startup(self):
        """Start a daophot session and attaches the inputImagePath's image.
        """
        # We start daophot from the working directory (the directory of the
        # input image.) All output will be placed in this directory. From
        # the user's perspective, the returned paths will still be relative
        # to the pipeline's base directory.
        startupCommand = '/bin/tcsh -c "cd %s;daophot"' % self.workDir
        self.daophot = pexpect.spawn(startupCommand)
        self.daophot.logfile = sys.stdout  # DEBUG
        self.daophot.expect("Command:")
        # print self.daophot.before
        self.set_option("WA", "-2")  # turn off extraneous printing
        self.attach('input_image')
    
    def shutdown(self):
        """Shutdown the daophot process."""
        self.daophot.sendline("exit")
        self.daophot = None
    
    def set_option(self, name, value):
        """Set the named option in daophot to a given value."""
        self.daophot.sendline("OPTION")
        self.daophot.expect(":")  # asks for the file with parameter values
        self.daophot.sendline("")  # accept the defaults
        self.daophot.expect("OPT>")
        self.daophot.sendline("%s=%s" % (name, value))
        self.daophot.expect("OPT>")
        self.daophot.sendline("")
        self.daophot.expect("Command:")
        print self.daophot.before
    
    def attach(self, image):
        """Attaches the given image to daophot. *image* will be resolved
        either as a name in the imageCache, or as a path. (Runs daophot
        *ATTACH*)
        
        By default, the attached image will be the last one attached (or the
        inputImagePath on the first run). But if *image* is specified, then
        it will be resolved in two steps
        1) If a name in the imageCache, that path will be used
        2) If not in the imageCache, then it will be used as a path itself
        
        FIXME change interface to work with both names and paths
        """
        imagePath = self._resolve_path(image, 'fits')
        self._set_last_path(imagePath, 'fits')
        
        command = "ATTACH %s" % imagePath
        self.daophot.sendline(command)
        self.daophot.expect("Command:")
    
    def find(self, nAvg=1, nSum=1, cooName=None, cooPath=None):
        """Runs the *FIND* command on the previously attached image.
        
        :param cooName: Set to have the coordinate path cached in the
            *self.cooCache*
        :param cooPath: Set as the filepath for the output coordinate file,
            otherwise a default path is made up.
        """
        cooPath = self._make_output_path(cooPath, cooName, "coo")
        self._name_path(cooName, cooPath, 'coo')
        self._set_last_path(cooPath, 'coo')
        
        self.daophot.sendline("FIND")
        # asks 'Number of frames averaged, summed:'
        self.daophot.expect(":")
        self.daophot.sendline("%i,%i" % (nAvg, nSum))
        # asks 'File for positions (default ???.coo):'
        self.daophot.expect(":")
        self.daophot.sendline(cooPath)
        self.daophot.expect("Are you happy with this?", timeout=60 * 20)
        # print self.daophot.before
        self.daophot.sendline("Y")
        self.daophot.expect("Command:")
    
    def apphot(self, coordinates, apRadPath=None, photOutputPath=None,
            photOutputName=None, options=None):
        """Run aperture photometry routine *PHOTOMETRY* in daophot.
        
        :param coordinates: refers to coordinates of stars from the find
            method; it is a string to be resolved either into a name in the
            cooCache, a filepath itself.
        :param apRadPath: is path to the aperture radii options file. This file
            must be in the working diretory. If None, then the default of
            'photo.opt' is assumed.
        """
        self.daophot.sendline("PHOTOMETRY")
        # asks for 'File with aperture radii (default photo.opt)'
        self.daophot.expect(":")
        
        if apRadPath is not None:
            self.daophot.sendline(os.path.basename(apRadPath))
        else:
            self.daophot.sendline("")  # assume default photo.opt file
        
        self.daophot.expect("PHO>")
        if options is not None:
            for optionName, optionValue in options.iteritems():
                self.daophot.sendline(optionName + "=" + optionValue)
                self.daophot.expect("PHO>")
        self.daophot.sendline("")
        
        # asks 'Input position file (default source/sky28k.coo):'
        self.daophot.expect(":")
        
        cooPath = self._resolve_path(coordinates, 'coo')
        self.daophot.sendline(cooPath)
        # asks 'Output file (default source/sky28k.ap):'
        self.daophot.expect(":")
        
        photOutputPath = self._make_output_path(photOutputPath,
                photOutputName, "ap")
        self._name_path(photOutputName, photOutputPath, 'ap')
        self._set_last_path(photOutputPath, 'ap')
        
        self.daophot.sendline(photOutputPath)
        self.daophot.expect("Command:", timeout=60 * 20)
    
    def pick_psf_stars(self, nStars, apPhot, starListPath=None,
            starListName=None, magLimit=99):
        """Picks *nStars* number of stars from the aperture photometry list
        that will be used as prototypes for making a PSF model;
        runs daophot *PICK*.
        
        :param apPhot: points to the aperture photometry list (made by
            apphot()). It is resolved into a name in apCache or a filepath to
            the .ap file
        :param nStars: is the number of stars to select, can be a str or int.
        :param starListPath: and starListName and the path/name that may
            be specified for the .lst file that lists the prototype psf stars.
        :param magLimit: is the limiting instrumental magnitude that can be
            used as a PSF prototype. Can be a str object.
        """
        magLimit = str(magLimit)
        nStars = str(int(nStars))
        apPhotPath = self._resolve_path(apPhot, 'ap')
        
        starListPath = self._make_output_path(starListPath,
                starListName, 'lst')
        self._name_path(starListName, starListPath, 'lst')
        self._set_last_path(starListPath, 'lst')
        
        self.daophot.sendline("PICK")
        # ask for input file name to .ap file
        self.daophot.expect(":")
        self.daophot.sendline(apPhotPath)
        # asks for 'Desired number of stars, faintest magnitude:'
        self.daophot.expect(":")
        self.daophot.sendline(",".join((nStars, magLimit)))
        # asks for output file path, .lst
        self.daophot.expect(":")
        # TODO implement output filepath
        self.daophot.sendline("")
        self.daophot.expect("Command:", timeout=60 * 10)
    
    def make_psf(self, apPhot, starList, psfPath=None, psfName=None):
        """Computes a PSF model with the daophot *PSF* command.
        
        :param apPhot: points to the aperture photometry list (made by
            apphot()). It is resolved into a name in apCache or a filepath
            to the .ap file
        :param starList: points to the psf prototype star list.
        
        :return: text output of fitting routine, path to the psf file and path
        to the neighbours file
        """
        apPhotPath = self._resolve_path(apPhot, 'ap')
        
        starListPath = self._resolve_path(starList, 'lst')
        
        psfPath = self._make_output_path(psfPath, psfName, 'psf')
        self._name_path(psfName, psfPath, 'psf')
        self._set_last_path(psfPath, 'psf')
        
        # make with the neighbours file name (.nei).
        # It always follows this form:
        fileRoot = os.path.splitext(psfPath)[0]
        neiPath = ".".join((fileRoot, 'nei'))
        self._set_last_path(neiPath, 'nei')
        if os.path.exists(neiPath):
            os.remove(neiPath)
        
        self.daophot.sendline("PSF")
        # asks for file with aperture phot results
        self.daophot.expect(":")
        self.daophot.sendline(apPhotPath)
        # asks for file with PSF prototype star list
        self.daophot.expect(":")
        self.daophot.sendline(starListPath)
        # asks for file for the psf output file
        self.daophot.expect(":")
        self.daophot.sendline(psfPath)
        # funny hack; pexpect has trouble here, but works
        # self.daophot.expect(".nei", timeout=120)
        
        # send a CR to make sure we're clean before leaving
        # self.daophot.sendline("")
        result = self.daophot.expect(["nei", "Failed to converge.",
            "Command:"], timeout=60 * 10)
        # save daophot's output of fit quality
        fittingText = self.daophot.before
        if result == 1 or result == 2:
            # failed to converge
            print "didn't converge. now what?"
            # raise PSFNotConverged
            return None, None, None
        
        # otherwise we should have good convergence
        print result,
        print "Ok convergence?"
        self.daophot.sendline("")
        self.daophot.expect("Command:")
        
        return fittingText, os.path.join(self.workDir, psfPath), \
            os.path.join(self.workDir, neiPath)
    
    def substar(self, substarList, psf, outputPath, keepers=None):
        """Subtracts stars in `substarList` from the attached image using the
        `psf` model.
        
        :param substarList: is a **path* to a photometry file of all stars
            that should be subtracted out of the image.
        :param psf: is a name/path resolved into a path to a PSF model.
        :param outputPath: is a **path** where the star-subtracted FITS image
            will be placed. Any existing file will be deleted.
        :param keepers: is a **path** (not a resolvable file name) to a
            listing stars that should be kept in the subtracted image.
            If `None`, then no stars are kept.
            
        :return: outputPath, relative to the pipeline.
        """
        psfPath = self._resolve_path(psf, 'psf')
        if os.path.exists(outputPath):
            os.remove(outputPath)
        
        self.daophot.sendline("SUBSTAR")
        self.daophot.expect(":")  # File with the PSF (*)
        self.daophot.sendline(os.path.basename(psfPath))
        self.daophot.expect(":")  # File with photometry (*)
        self.daophot.sendline(os.path.basename(substarList))
        # print "send substarList"
        # print self.daophot.interact()
        self.daophot.expect("in\?")  # Do you have stars to leave in
        if keepers is not None:
            self.daophot.sendline("Y")
            # print self.daophot.before
            self.daophot.expect(":")  # File with star list (*)
            self.daophot.sendline(os.path.basename(keepers))
        else:
            self.daophot.sendline("N")
        self.daophot.expect(":")  # Name for subtracted image (*)
        self.daophot.sendline(os.path.basename(outputPath))
        self.daophot.expect("Command:", timeout=60 * 10)
        
        return outputPath
    
    def get_path(self, name, ext):
        """Returns the named path of type ext. The path will be relative
        to the pipeline's base... as the user would expect."""
        return os.path.join(self.workDir, self._resolve_path(name, ext))
    
    def _resolve_path(self, path, ext):
        """Resolves path into a path to the given type (via ext extension)
        of file if it is name. Or if it is a path already, that
        path will be passed through. The returned path is relative to the
        workDir (working directory) of this Daophot.
        """
        print path,
        print ext
        try:
            resolvedPath = self.pathCache[ext][path]
        except:
            print "This is a path"
            print path
            resolvedPath = os.path.basename(path)
        return resolvedPath
    
    def _name_path(self, name, path, ext):
        """Adds the path of type ext(ention) to its cache under given name,
        if the name is not None.
        """
        if name is not None:
            self.pathCache[ext][name] = path
    
    def _set_last_path(self, path, ext):
        """Makes the path be filed under 'last' in its type's cache."""
        self.pathCache[ext]['last'] = path
    
    def _make_output_path(self, path, name, ext):
        """Forms an output file path. If path is None, then a path is made
        using the name. If both path and name are None, then a path is formed
        from the inputImagePath and the filename extension *ext*.
        
        The path is force to be relative to `workDir`.
        """
        if path is None:
            # make a default ap photometry output file path
            fileRoot = os.path.splitext(
                    os.path.basename(self.inputImagePath))[0]
            if name is not None:
                fileRoot = "_".join((fileRoot, name))
            path = ".".join((fileRoot, ext))
        else:
            path = os.path.basename(path)
        
        fullpath = os.path.join(self.workDir, path)
        if os.path.exists(fullpath):
            print "removing existing %s" % fullpath
            os.remove(fullpath)
        
        return path
