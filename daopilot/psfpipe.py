import os
import shutil
import re
import glob
import numpy
import pyfits
from astLib import astWCS

import owl.region
import owl.twomicron
import owl.Match

from daophot import Daophot
from allstar import Allstar


class PSFFactory(object):
    """Factory class for creating PSFs from a single image.
    """
    def __init__(self, workDir, apRadPath):
        super(PSFFactory, self).__init__()
        self.workDir = workDir
        if not os.path.exists(self.workDir): os.makedirs(workDir)

        # Install aperture radius config file
        self.apRadPath = os.path.join(self.workDir, "photo.opt")
        if apRadPath != self.apRadPath:
            if os.path.exists(self.apRadPath): os.remove(apRadPath)
            shutil.copy(apRadPath, self.apRadPath)
    
    def make(self, imageName, imagePath, flagPath, band, maxVarPSF,
            runAllstar=False, findHiddenStars=False, clean=False):
        """Makes the PSF model.
        
        :param maxVarPSF: the maximum degrees of freedom in the PSF. Maximum
            is 2.
        :param runAllstar: set to True if you want an allstar star-subtracted
            image documenting each step of the psf subtraction process.
        """
        self.imageName = imageName
        self.imagePath = imagePath
        self.flagPath = flagPath
        self.band = band
        
        self.findHiddenStars = findHiddenStars
        
        self.daophot = Daophot(self.imagePath)
        
        # Initial DAOPHOT run
        # FIND
        self.daophot.set_option('VA', '-1')  # full analytic psf
        self.daophot.find(nAvg=1, nSum=1)
        coordFilePath = self.daophot.get_path('last', 'coo')
        findN, findX, findY = owl.dao.parseCoordFile(coordFilePath)
        findPoints = owl.region.PointList()
        findPoints.setFrame('image')
        findPoints.setPoints(findX, findY, size=6, shapes="circle",
                labels=None, colours="yellow")
        findPoints.writeTo(os.path.join(self.workDir,
                self.imageName + "_find.reg"))
        
        # PICK PSF stars
        self.daophot.apphot(coordinates='last',
                apRadPath=os.path.basename(self.apRadPath))
        apFilePath = self.daophot.get_path('last', 'ap')
        self.daophot.pickPSFStars(100, apPhot='last')
        
        psfN, psfX, psfY = owl.dao.parseCoordFile(
            self.daophot.get_path('last', 'lst'))
        psfPoints = owl.region.PointList()
        psfPoints.setFrame('image')
        psfPoints.setPoints(psfX, psfY, size=15, shapes="diamond",
                labels=None, colours="red")
        psfPoints.writeTo(os.path.join(self.workDir,
            self.imageName + "_psf.reg"))
        
        # Make custom picks
        picker = StarPicker(self.daophot, 'last', self.imagePath)
        picker.useDaophotPicks()
        if self.flagPath is not None:
            picker.filterOnFlagMap(self.flagPath)
        picker.filterBright2MASSByDistance(40., 14., self.band)
        pickPath = os.path.join(self.workDir, self.imageName + "rev.lst")
        picker.write(pickPath)
        picker.writeRegions(os.path.join(self.workDir,
            self.imageName + "_psfrev.reg"))
        
        # Make PSF, run allstar
        fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                starList=pickPath, psfName='init')
        alsPath, alsStarSubPath, neiSubPath = self._makeAllstarPaths("init")
        if runAllstar:
            allstar = Allstar(self.imagePath,
                self.daophot.get_path('last', 'psf'),
                self.daophot.get_path('last', 'ap'), alsPath, alsStarSubPath)
            allstar.run()
        
        # iterative DAOPHOT runs with increasing psf variability
        # psfStarListPath = self.daophot.get_path('last', 'lst')
        for varPSF in range(0, maxVarPSF + 1):
            itername = "var%i" % varPSF
            try:
                self._iteratePSF(varPSF, picker, neiPath, runAllstar,
                        name=itername)
            except PSFNotConverged:
                varPSF = -1
                self._makeAnalyticPSF(picker)
        
        # make final psf on clean image, keeping last-used varPSF
        success = self._iteratePSF(varPSF, picker, neiPath, runAllstar,
                name='fin')
        print success
        
        # Get path to the final psf
        psfPath = self.daophot.get_path("fin", "psf")
        
        self.daophot.shutdown()
        
        if clean:
            self._clean()
        
        return (psfPath, pickPath, coordFilePath, apFilePath)
    
    def _makeAnalyticPSF(self, picker):
        """This is a bailout method to return the path to the analytic PSF.
        This is called whenever the empirical PSFs fail to converge."""
        print "FALLING BACK TO ANALYTIC PSF"
        self.daophot.set_option("VA", "-1")
        pickPath = picker.getOutputPath()
        fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                starList=pickPath, psfName='bail')
        self.daophot.shutdown()
        return psfPath
    
    def _makeAllstarPaths(self, itername):
        """Makes paths for the Allstar photometry file (.als) and the
        star-subtracted FITS file automatically based on the input file path.
        """
        imageRoot = os.path.splitext(self.imagePath)[0]
        alsPath = ".".join((imageRoot, itername, "als"))
        alsStarSubPath = "_".join((imageRoot, itername, "als.fits"))
        neiSubPath = "_".join((imageRoot, itername, "subnei.fits"))
        return alsPath, alsStarSubPath, neiSubPath
    
    def _iteratePSF(self, varPSF, starPicker, neiPath, runAllstar, name=None):
        """Performs a recipe of
        * subtract neighbouring stars
        * fit PSF
        * make new allstar catalog
        """
        if name is None:
            name = str(varPSF)
        alsPath, alsStarSubPath, neiSubPath = self._makeAllstarPaths(name)
        pickPath = starPicker.getOutputPath()

        # had PSF fit be repeated; will reset to false if no stars are culled
        repeat = True
        while repeat:
            neiSubPath = self.daophot.substar(neiPath, 'last', neiSubPath,
                    keepers=pickPath)
            self.daophot.set_option('VA', int(varPSF))
            self.daophot.attach(neiSubPath)  # use the nei-subtracted image
            fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                    starList=pickPath, psfName=name)
            if starPicker.cullWithFitResults(fitText):
                repeat = True
                print "repeating PSF fit after culling stars"
            else:
                repeat = False
        
        if runAllstar:
            allstar = owl.dao.Allstar(self.imagePath,
                    self.daophot.get_path('last', 'psf'),
                    self.daophot.get_path('last', 'ap'),
                    alsPath, alsStarSubPath)
            allstar.run()
        
        if self.findHiddenStars:
            self.detectHiddenStars(self.daophot.get_path('last', 'psf'),
                    self.daophot.get_path('last', 'ap'),
                    alsPath, alsStarSubPath)
    
    def detectHiddenStars(self, psfPath, apPhotPath, alsPath, alsStarSubPath):
        """Runs allstar with the most current psf model; runs daophot find
        on that star-subtracted image and attempts to uncover new stars.
        """
        allstar = Allstar(self.imagePath, psfPath, apPhotPath,
                alsPath, alsStarSubPath)
        allstar.run()
        
        starSubDaophot = owl.dao.Daophot(alsStarSubPath)
        starSubDaophot.startup()
        starSubDaophot.find()
        starSubDaophot.apphot('last', apRadPath="wirphoto.opt")
        newApPath = starSubDaophot.get_path('last', 'ap')
        
        originalApCatalog = owl.dao.ApPhotCatalog2()
        originalApCatalog.open(apPhotPath)
        
        newApCatalog = owl.dao.ApPhotCatalog2()
        newApCatalog.open(newApPath)
        
        imageRoot = os.path.splitext(alsStarSubPath)[0]
        regPath = imageRoot + "_hidden.reg"
        newApCatalog.writeImageRegions(regPath)
        print "==== Detected %s hidden stars ====" % newApCatalog.nStars
        
        originalApCatalog.appendCatalog(newApCatalog)
        originalApCatalog.write(apPhotPath)  # write new catalog in place!
        
        starSubDaophot.shutdown()
    
    def _clean(self):
        """Uses a simple recipe to guess/try the names of files that should be
        deleted to save space.
        """
        imageRoot = os.path.splitext(self.imagePath)[0]
        globPaths = glob.glob(imageRoot + "*.fits")
        
        for path in globPaths:
            if (path == self.imagePath) | (path == self.flagPath):
                print path
                continue
            else:
                os.remove(path)
        
        for path in glob.glob(imageRoot + "_var0*"):
            os.remove(path)
        for path in glob.glob(imageRoot + "_var1*"):
            os.remove(path)
        for path in glob.glob(imageRoot + "_var2*"):
            os.remove(path)
        for path in glob.glob(imageRoot + "*.nei"):
            os.remove(path)
        for path in glob.glob(imageRoot + "_init*"):
            os.remove(path)
        for path in glob.glob(imageRoot + "*.als"):
            os.remove(path)
        # for path in glob.glob(imageRoot+"_.lst"):
            # os.remove(path)
        # for path in glob.glob(imageRoot+"_.coo"):
            # os.remove(path)
        for path in glob.glob(imageRoot + "_fin_als.*"):
            os.remove(path)
        for path in glob.glob(imageRoot + "*.reg"):
            os.remove(path)
        os.remove(imageRoot + ".lst")
        # os.remove(imageRoot+".coo")


class PSFNotConverged(Exception): pass


class StarPicker(object):
    """StarPicker is intended as a replacement for the built-in DAOPHOT/PICK.
    """
    def __init__(self, daophot, daophotName, inputImagePath):
        """
        :param daophot: a `Daophot` instance of the image being worked on.
        :param daophotName: name that has the cached results (from, e.g. FIND,
            PHOTOMETRY) to be used to build the star lists. All types of
            results should be cached in `daophot` under the same name.
        """
        super(StarPicker, self).__init__()
        self.daophot = daophot
        self.daophotName = daophotName
        self.inputImagePath = inputImagePath
        self.psc = None  # 2MASS point source catalog
        self.outputPath = None  # where the lst will be saved
        
        # Load the aperture photometry of the star list
        self.apCatalog = owl.dao.ApPhotCatalog(daophot, daophotName)
        # By default, accept all stars as PSF candidates
        self.candidates = self.apCatalog.getStarIDs()
        print "There are %i candidates on init" % len(self.candidates)
    
    def getOutputPath(self):
        return self.outputPath
    
    def write(self, outputPath):
        print "There are %i candidates on write" % len(self.candidates)
        """Writes the list of selected PSF model stars to disk."""
        x = [self.apCatalog.stars[star]['x'] for star in self.candidates]
        y = [self.apCatalog.stars[star]['y'] for star in self.candidates]
        mag = [self.apCatalog.stars[star]['mag'] for star in self.candidates]
        magErr = [self.apCatalog.stars[star]['mag_err']
                for star in self.candidates]
        pickCatalog = owl.dao.PickCatalog()
        pickCatalog.setStars(self.candidates, x, y, mag, magErr)
        pickCatalog.write(outputPath)
        self.outputPath = outputPath
    
    def writeRegions(self, outputPath):
        """Writes a DS9-compatible .reg file with the selected PSF model stars.
        """
        x = [self.apCatalog.stars[star]['x'] for star in self.candidates]
        y = [self.apCatalog.stars[star]['y'] for star in self.candidates]
        names = [str(star) for star in self.candidates]
        
        psfPoints = owl.region.PointList()
        psfPoints.setFrame('image')
        psfPoints.setPoints(x, y, size=15, shapes="circle", labels=names,
                colours="cyan")
        psfPoints.writeTo(outputPath)
    
    def useDaophotPicks(self):
        """Whittles down the candidate list to just those selected by DAOPHOT
        PICK.
        """
        pickCatalog = owl.dao.PickCatalog()
        pickCatalog.readFromDaophot(self.daophot, self.daophotName)
        self.candidates = pickCatalog.getStarIDs()
        print "There are %i candidates on useDaophot" % len(self.candidates)
    
    def filterOnFlagMap(self, flagPath):
        """Applies the flagmap to filtering the PSF template stars. Any star
        whose centroid lies upon a flagged (>0) pixel will be rejected
        from candidacy.
        
        :param flagPath: is the **filepath** to the flag image (not the
            flagName!; this is done because DAOPHOT works on single extension
            images; I don't have a good way of referring to a certain extension
            of a flag image yet.)
        """
        flagFITS = pyfits.open(flagPath)
        for star in self.candidates:
            x = int(self.apCatalog.stars[star]['x'])
            y = int(self.apCatalog.stars[star]['y'])
            if flagFITS[0].data[y, x] > 0:
                self.candidates.remove(star)
        flagFITS.close()
        print "There are %i candidates on flagmap filter" % \
                len(self.candidates)
    
    def _get2MASS(self):
        """Sets the self.psc catalog with 2MASS stars in the input image frame.
        """
        inputFITS = pyfits.open(self.inputImagePath)
        header = inputFITS[0].header
        wcs = astWCS.WCS(header, mode='pyfits')
        alphaNW, deltaNW = wcs.pix2wcs(2048, 2048)
        alphaSW, deltaSW = wcs.pix2wcs(2048, 1)
        alphaSE, deltaSE = wcs.pix2wcs(1, 1)
        alphaNE, deltaNE = wcs.pix2wcs(1, 2048)
        raMin = min(alphaNW, alphaSW, alphaNE, alphaNW)
        raMax = max(alphaNW, alphaSW, alphaNE, alphaNW)
        decMin = min(deltaNW, deltaSW, deltaNE, deltaNW)
        decMax = max(deltaNW, deltaSW, deltaNE, deltaNW)
        catalog2MASS = owl.twomicron.Catalog2MASS()
        self.psc = catalog2MASS.getStarsInArea(raMin, raMax, decMin, decMax)
        inputFITS.close()
    
    def filterBright2MASSByDistance(self, radius, magThreshold, band):
        """Rejects stars that are a certain distance from bright 2MASS stars
        that are known by the user to be saturated. This might be the
        candidate star itself.
        
        I propose that 13th magnitude is a good cut-off for a star that you'd
        want to be close to. Need to check this...
        
        :param radius: exclusion zone of bright 2MASS stars around the
            PSF candidate, in pixels.
        :param band: can either be "J" or "Ks"
        """
        radius = radius ** 2.  # work in squared distances
        
        if band == "J":
            magKey = "jmag"
        elif band == "Ks":
            magKey = "kmag"
        else:
            magKey = None
        
        brightRA = []
        brightDec = []
        
        # Get 2MASS stars and magnitudes in the native frame of the input image
        if self.psc is None:
            self._get2MASS()
        
        # Find stars brighter than threshold
        n2MASS = len(self.psc['ra'])
        for i in xrange(n2MASS):
            if self.psc[magKey][i] < magThreshold:
                # convert this to image frame with pointToImageFrame...
                brightRA.append(self.psc['ra'][i])
                brightDec.append(self.psc['dec'][i])
        nBrights = len(brightRA)
        brightX = numpy.zeros(nBrights)
        brightY = numpy.zeros(nBrights)
        inputFITS = pyfits.open(self.inputImagePath)
        header = inputFITS[0].header
        wcs = astWCS.WCS(header, mode='pyfits')
        for i in xrange(nBrights):
            x, y = wcs.wcs2pix(brightRA[i], brightDec[i])
            brightX[i] = x
            brightY[i] = y
        
        # Ask candidate if it is within radius pixels of a bright star
        for star in self.candidates:
            x = int(self.apCatalog.stars[star]['x'])
            y = int(self.apCatalog.stars[star]['y'])
            dist = (x - brightX) ** 2. + (y - brightY) ** 2.
            # print "The closest bright star is %.2f pixels" % dist.min()
            if len(numpy.where(dist < radius)[0]) > 0:
                # there are bright stars nearby; need to delete this candidate
                # print "There is a bright 2MASS star on/near %i" % star
                # print self.candidates.index(star)
                self.candidates.remove(star)
        
        print "There are %i candidates on 2MASS filter" % len(self.candidates)
    
    def filterByNeighbours(self):
        """Rejects stars that have neighbours within a certain distance.
        
        This method is different from `filterBright2MASSByDistance` as that one
        rejects a star if any bright 2MASS star exists within a radius about
        a candidate---it may be the candidate itself. This method simply asks
        if this given star has a bright neighbour and doesn't
        *self obliterate*.
        """
        pass
    
    def cullWithFitResults(self, fitText):
        """Removes stars flagged in `fitText`, which is returned by `daophot`
        during the psf fitting process. The culled stars are removed from
        `self.candidates`, and a new .lst file is automatically save to
        `self.outputPath`, where the previous .lst had been saved. This method
        is a way of using daophot feedback to improve the PSF star list.
        
        There are two types of flags that are detected. First:
        *. '   2182 is not a good star.' will cause star 2182 to be removed
        *. '   2051  0.054      1648  0.059      1342  0.068       471  0.061      1522  0.121 ?'
            will cause star 1522 to be removed.
        
        :return: `True` if stars were culled, `False` otherwise.
        """
        
        badStars = []  # list of star IDs to be removed
        lines = fitText.split("\n")
        
        # filter for e.g. " 2182 is not a good star."
        for line in lines:
            matches = re.match("[ \n\t]*([0-9]*) is not a good star.[\n]*",
                    line)
            if matches:
                for group in matches.groups():
                    badStars.append(int(group))
        
        # filter for the "?" flags
        # character positions for the '?' flags
        flagPositions = (15, 32, 49, 66, 83)
        for line in lines:
            print line
            try:
                for flagPos in flagPositions:
                    if (line[flagPos] == "?") or (line[flagPos] == "*"):
                        # print "bad on %i" % flagPos
                        badStars.append(int(line[flagPos - 8 - 7:flagPos - 8]))
            except:
                continue
        
        print "bad stars:",
        print badStars
        
        # Cull and save new star list
        for badStar in badStars:
            self.candidates.remove(badStar)
        self.write(self.outputPath)
        
        if len(badStars) > 0:
            return True
        else:
            return False
