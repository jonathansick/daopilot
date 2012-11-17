import os
import shutil
import re
import glob
import pyfits

import regionio
import catalogio

from daophot import Daophot
from allstar import Allstar


class PSFFactory(object):
    """Factory class for creating PSFs from a single image.

    :param picker: Instance of PSFPicker, or subclass. Used to refine the
        PSF model star list.
    """
    def __init__(self, workDir, apRadPath, picker):
        super(PSFFactory, self).__init__()
        self.workDir = workDir
        if not os.path.exists(self.workDir): os.makedirs(workDir)
        self.picker = picker

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
        cooCatalog = catalogio.CoordCatalog()
        cooCatalog.open(coordFilePath)
        findX = cooCatalog.stars['x']
        findY = cooCatalog.stars['y']
        findPoints = regionio.PointList()
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
        
        # Make custom picks for PSF models
        # TODO force use of lst; or allow the picker to work directly
        # from apphot?
        pickPath = os.path.join(self.workDir, self.imageName + "_rev.lst")
        pickRegPath = os.path.join(self.workDir,
                self.imageName + "_psfrev.reg")
        self.picker(self.daophot, 'last', self.imagePath, pickPath,
                regionPath=pickRegPath)
        
        # Make PSF, run allstar
        fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                starList=pickPath, psfName='init')
        alsPath, alsStarSubPath, neiSubPath = self._make_allstar_paths("init")
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
                self._iterate_psf(varPSF, neiPath, runAllstar,
                        name=itername)
            except PSFNotConverged:
                varPSF = -1
                self._make_analytic_psf()
        
        # make final psf on clean image, keeping last-used varPSF
        success = self._iterate_psf(varPSF, neiPath, runAllstar,
                name='fin')
        print success
        
        # Get path to the final psf
        psfPath = self.daophot.get_path("fin", "psf")
        
        self.daophot.shutdown()
        
        if clean:
            self._clean()
        
        return (psfPath, pickPath, coordFilePath, apFilePath)
    
    def _make_analytic_psf(self):
        """This is a bailout method to return the path to the analytic PSF.
        This is called whenever the empirical PSFs fail to converge."""
        print "FALLING BACK TO ANALYTIC PSF"
        self.daophot.set_option("VA", "-1")
        pickPath = self.picker.get_output_path()
        fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                starList=pickPath, psfName='bail')
        self.daophot.shutdown()
        return psfPath
    
    def _make_allstar_paths(self, itername):
        """Makes paths for the Allstar photometry file (.als) and the
        star-subtracted FITS file automatically based on the input file path.
        """
        imageRoot = os.path.splitext(self.imagePath)[0]
        alsPath = ".".join((imageRoot, itername, "als"))
        alsStarSubPath = "_".join((imageRoot, itername, "als.fits"))
        neiSubPath = "_".join((imageRoot, itername, "subnei.fits"))
        return alsPath, alsStarSubPath, neiSubPath
    
    def _iterate_psf(self, varPSF, neiPath, runAllstar, name=None):
        """Performs a recipe of
        * subtract neighbouring stars
        * fit PSF
        * make new allstar catalog
        """
        if name is None:
            name = str(varPSF)
        alsPath, alsStarSubPath, neiSubPath = self._make_allstar_paths(name)
        pickPath = self.picker.get_output_path()

        # had PSF fit be repeated; will reset to false if no stars are culled
        repeat = True
        while repeat:
            neiSubPath = self.daophot.substar(neiPath, 'last', neiSubPath,
                    keepers=pickPath)
            self.daophot.set_option('VA', int(varPSF))
            self.daophot.attach(neiSubPath)  # use the nei-subtracted image
            fitText, psfPath, neiPath = self.daophot.make_psf(apPhot='last',
                    starList=pickPath, psfName=name)
            if self.picker.cull_with_fit_results(fitText):
                repeat = True
                print "repeating PSF fit after culling stars"
            else:
                repeat = False
        
        if runAllstar:
            allstar = Allstar(self.imagePath,
                    self.daophot.get_path('last', 'psf'),
                    self.daophot.get_path('last', 'ap'),
                    alsPath, alsStarSubPath)
            allstar.run()
        
        if self.findHiddenStars:
            self._detect_hidden_stars(self.daophot.get_path('last', 'psf'),
                    self.daophot.get_path('last', 'ap'),
                    alsPath, alsStarSubPath)
    
    def _detect_hidden_stars(self, psfPath, apPhotPath, alsPath,
            alsStarSubPath):
        """Runs allstar with the most current psf model; runs daophot find
        on that star-subtracted image and attempts to uncover new stars.
        """
        allstar = Allstar(self.imagePath, psfPath, apPhotPath,
                alsPath, alsStarSubPath)
        allstar.run()
        
        starSubDaophot = Daophot(alsStarSubPath)
        starSubDaophot.startup()
        starSubDaophot.find()
        starSubDaophot.apphot('last', apRadPath="wirphoto.opt")
        newApPath = starSubDaophot.get_path('last', 'ap')
        
        originalApCatalog = catalogio.ApPhotCatalog2()
        originalApCatalog.open(apPhotPath)
        
        newApCatalog = catalogio.ApPhotCatalog2()
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
    def __init__(self):
        super(StarPicker, self).__init__()

    def __call__(self, daophot, daophotName, outputPath, regionPath=None):
        """Run the PSF star picking pipeline.
        
        :param daophot: a `Daophot` instance of the image being worked on.
        :param daophotName: name that has the cached results (from, e.g. FIND,
            PHOTOMETRY) to be used to build the star lists. All types of
            results sould be cached in `daophot` under the same name.
        """
        self.daophot = daophot
        self.daophotName = daophotName
        self.outputPath = None  # where the lst will be saved
        
        # Load the aperture photometry of the star list
        self.apCatalog = catalogio.ApPhotCatalog(daophot, daophotName)
        # By default, accept all stars as PSF candidates
        self.candidates = self.apCatalog.getStarIDs()
        print "There are %i candidates on init" % len(self.candidates)

        self.user_processing()

        self.write(outputPath)

        if regionPath is not None:
            self.write_regions(regionPath)

    def user_processing(self):
        """This method should be implemented by the user."""
        pass
    
    def get_output_path(self):
        return self.outputPath
    
    def write(self, outputPath):
        print "There are %i candidates on write" % len(self.candidates)
        """Writes the list of selected PSF model stars to disk."""
        x = [self.apCatalog.stars[star]['x'] for star in self.candidates]
        y = [self.apCatalog.stars[star]['y'] for star in self.candidates]
        mag = [self.apCatalog.stars[star]['mag'] for star in self.candidates]
        magErr = [self.apCatalog.stars[star]['mag_err']
                for star in self.candidates]
        pickCatalog = catalogio.PickCatalog()
        pickCatalog.setStars(self.candidates, x, y, mag, magErr)
        pickCatalog.write(outputPath)
        self.outputPath = outputPath
    
    def write_regions(self, outputPath):
        """Writes a DS9-compatible .reg file with the selected PSF model stars.
        """
        x = [self.apCatalog.stars[star]['x'] for star in self.candidates]
        y = [self.apCatalog.stars[star]['y'] for star in self.candidates]
        names = [str(star) for star in self.candidates]
        
        psfPoints = regionio.PointList()
        psfPoints.setFrame('image')
        psfPoints.setPoints(x, y, size=15, shapes="circle", labels=names,
                colours="cyan")
        psfPoints.writeTo(outputPath)
    
    def filter_on_flagmap(self, flagPath):
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
    
    def cull_with_fit_results(self, fitText):
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
