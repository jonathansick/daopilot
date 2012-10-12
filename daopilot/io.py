#!/usr/bin/env python
# encoding: utf-8
"""
Classes for reading and writing Daophot files.

2012-05-05 - Created by Jonathan Sick
"""

import os
import numpy as np

import region  # TODO change this to pyregion/astropy package


class DaoCatalogBase(object):
    """Base class for the suite of DAOPHOT I/O catalogs."""
    def __init__(self):
        super(DaoCatalogBase, self).__init__()
        self.stars = None
        self.headerText = None
        self.nStars = 0
        self.nHeaderLines = 2

    def open(self, path):
        """docstring for open"""
        catfile = open(path)
        headerLines, dataLines = self._split_header(catfile)
        catfile.close()
        self.parse(dataLines)
        self.headerText = "".join(headerLines)
    
    def _split_header(self, f):
        """Given a catalog file descriptor, returns lists of text lines, split
        between the header and data components.
        """
        lines = f.readlines()
        return lines[0:self.nHeaderLines + 1], lines[self.nHeaderLines + 1:]
    
    def get_header(self):
        return self.headerText
    
    def set_header(self, headerText):
        self.headerText = headerText
    
    def write_regions(self, outputPath, markersize=10, markercolour='red',
            marker="circle"):
        """Writes a .reg file in the image frame with all sources
        
        .. todo:: Port this to pyregion/astropy package
        """
        pointList = region.PointList()
        pointList.setFrame('image')
        pointList.setPoints(self.stars['x'], self.stars['y'], shapes=marker,
                colours=markercolour, size=markersize)
        pointList.writeTo(outputPath)
    
    def append_catalog(self, newCatalog):
        """Appends a catalog to the end of the current catalog. The serial
        numbers of stars in the `newCatalog` are updated to be continuous
        with the current catalog's.
        """
        maxExistingID = self.stars['id'].max()
        newCatalog.stars['id'] = newCatalog.stars['id'] + maxExistingID
        self.stars = np.concatenate((self.stars, newCatalog.stars))
        self.nStars = len(self.stars)
    
    def write(self, outputPath):
        """Saves the catalog data in the .ap format to `outputPath`."""
        catalogLines = self.make_catalog_lines()
        catalogText = self.headerText + "\n".join(catalogLines)
        
        if os.path.exists(outputPath):
            os.remove(outputPath)
        
        f = open(outputPath, 'w')
        f.write(catalogText)
        f.close()
    
    def right_align_int(self, number, length):
        """docstring for integerLine"""
        s = "%i" % number
        if len(s) < length:
            s = " " * (length - len(s)) + s
        return s
    
    def right_align_F3(self, number, lengthBefore):
        """docstring for right_align_int3"""
        s = "%.3f" % number
        if len(s) < (lengthBefore + 4):
            s = " " * (lengthBefore + 4 - len(s)) + s
        return s
    
    def right_align_F2(self, number, lengthBefore):
        """docstring for right_align_int2"""
        s = "%.3f" % number
        if len(s) < (lengthBefore + 3):
            s = " " * (lengthBefore + 3 - len(s)) + s
        return s
    
    def mag_str(self, number):
        """docstring for mag_str"""
        if number >= 99.999:
            s = "99.999"
        else:
            s = "%.3f" % number
        
        if len(s) < 6:
            s = " " * (6 - len(s)) + s
        return s
    
    def mag_err_str(self, number):
        """docstring for mag_err_str"""
        if number >= 9.9999:
            s = "9.9999"
        else:
            s = "%.4f" % number
        
        if len(s) < 6:
            s = " " * (6 - len(s)) + s
        return s


class CoordCatalog(DaoCatalogBase):
    """For managing (reading/writing) .coo files, like produced by
    daophot FIND
    """
    def __init__(self):
        super(CoordCatalog, self).__init__()
        self.fullCatalog = False  # True if it has 6 data records per line
        self.nHeaderLines = 2
        self.dt = np.dtype([('id', np.uint), ('x', np.float32),
            ('y', np.float32), ('mag', np.float32), ('sharpness', np.float32),
            ('roundness', np.float32), ('marginal_roundness', np.float32)])
    
    def parse(self, dataLines):
        """docstring for parse"""
        self.nStars = len(dataLines)
        self.stars = np.empty(self.nStars, dtype=self.dt)
        
        for i, line in enumerate(dataLines):
            items = line.split()
            self.stars['id'][i] = int(items[0])
            self.stars['x'][i] = float(items[1])
            self.stars['y'][i] = float(items[2])
            if len(items) > 3:
                # expect a full catalog
                self.stars['mag'][i] = float(items[3])
                self.stars['sharpness'][i] = float(items[4])
                self.stars['roundness'][i] = float(items[5])
                self.stars['marginal_roundness'][i] = float(items[6])
                self.fullCatalog = True
    
    def set_stars(self, newIDs, newX, newY, newMag, newSharpness, newRoundness,
            newMarginalRoundness):
        """docstring for set_stars"""
        self.nStars = len(newIDs)
        self.stars = np.empty(self.nStars, dtype=self.dt)
        self.stars['id'] = newIDs
        self.stars['x'] = newX
        self.stars['y'] = newY
        if newMag is not None:
            self.stars['mag'] = newMag
        if newSharpness is not None:
            self.stars['sharpness'] = newSharpness
        if newRoundness is not None:
            self.stars['roundness'] = newRoundness
        if newMarginalRoundness is not None:
            self.stars['marginal_roundness'] = newMarginalRoundness
    
    def make_catalog_lines(self):
        catalogLines = []
        for i in xrange(self.nStars):
            line = "% *i % *.3f % *.3f %.3f %.3f %.3f %.3f" \
                    % (3, self.stars['id'][i], 4, self.stars['x'][i],
                    4, self.stars['y'][i], self.stars['mag'][i],
                    self.stars['sharpness'][i],
                    self.stars['roundness'][i],
                    self.stars['marginal_roundness'][i])
            catalogLines.append(line)
        return catalogLines


class ApPhotCatalog(DaoCatalogBase):
    """A revised class for reading .ap catalogs produced by Daophot:Photometry.
    """
    def __init__(self):
        super(ApPhotCatalog, self).__init__()
        self.nHeaderLines = 3
    
    def parse(self, dataLines):
        """docstring for parse"""
        dt = np.dtype([('id', np.uint), ('x', np.float32), ('y', np.float32),
            ('mag', np.float32), ("modal_sky", np.float32),
            ("sky_sigma", np.float32), ("sky_skew", np.float32),
            ("mag_err", np.float32)])
        
        # In .ap catalogs, each star has data on two lines; lets split text
        # into this grouping so that we can enumerate over the stars themselves
        dataText = "".join(dataLines)
        # double newline separates stars
        groupedStarLines = dataText.split("\n\n")
        
        self.nStars = len(groupedStarLines)
        self.stars = np.empty(self.nStars, dtype=dt)
        
        for i, groupedStarLine in enumerate(groupedStarLines):
            lines = groupedStarLine.split("\n")
            firstLineItems = lines[0].split()
            secondLineItems = lines[1].split()
            
            self.stars['id'][i] = int(firstLineItems[0])
            self.stars['x'][i] = float(firstLineItems[1])
            self.stars['y'][i] = float(firstLineItems[2])
            self.stars['mag'][i] = float(firstLineItems[3])
            
            self.stars['modal_sky'][i] = float(secondLineItems[0])
            self.stars['sky_sigma'][i] = float(secondLineItems[1])
            self.stars['sky_skew'][i] = float(secondLineItems[2])
            self.stars['mag_err'][i] = float(secondLineItems[3])
    
    def make_catalog_lines(self):
        catalogLines = []
        for i in xrange(self.nStars):
            firstLine = "%s %s %s %s" % \
                    (self.right_align_int(self.stars[i]['id'], 7),
                    self.right_align_F3(self.stars[i]['x'], 4),
                    self.right_align_F3(self.stars[i]['y'], 4),
                    self.mag_str(self.stars[i]['mag']))
            secondLine = "%s %s %s %s" % \
                    (self.right_align_F3(self.stars[i]['modal_sky'], 10),
                    self.right_align_F2(self.stars[i]['sky_sigma'], 2),
                    self.right_align_F2(self.stars[i]['sky_skew'], 2),
                    self.mag_err_str(self.stars[i]['mag_err']))
            catalogLines.append(firstLine + "\n" + secondLine + "\n")
        return catalogLines


class PickCatalog(object):
    """Reads the .lst catalogs produced by DAOPHOT's PICK routine."""
    def __init__(self):
        super(PickCatalog, self).__init__()
        self.stars = None
    
    def read_from_daophot(self, daophot, lstName):
        """Reads the named list from teh daophot instance."""
        lstPath = daophot.getPath(lstName, 'lst')
        self.read(lstPath)
    
    def read(self, lstPath):
        """Loads the .lst file at lstPath into the instance memory."""
        self.stars = {}
        
        f = open(lstPath, 'rU')
        for lineNumber, line in enumerate(f):
            if lineNumber < 3:
                continue  # I presume these are header lines
            items = line.split()
            serial = int(items[0])
            self.stars[serial] = {'x': float(items[1]),
                                  'y': float(items[2]),
                                  'mag': float(items[3]),
                                  'mag_err': float(items[4])}
        f.close()
    
    def get_star_ids(self):
        """Returns a list of the ID serial numbers of all stars in the catalog.
        """
        return self.stars.keys()
    
    def set_stars(self, serial, x, y, mag, magErr):
        """Sets the self.stars dictionary by providing lists of the star ID,
        position and magnitudes."""
        self.stars = {}
        for i, star in enumerate(serial):
            self.stars[star] = {'x': x[i], 'y': y[i], 'mag': mag[i],
                    'mag_err': mag[i]}
    
    def write(self, outputPath):
        """Writes a .lst file to the output path, based on the stars in the
        instance.
        """
        if os.path.exists(outputPath):
            os.remove(outputPath)
        f = open(outputPath, 'w')
        for idnum, star in self.stars.iteritems():
            line = "% 8i %.3f %.3f %.3f %.4f\n" % (idnum, star['x'], star['y'],
                    star['mag'], star['mag_err'])
            f.write(line)
        f.close()
    
    def write_regions(self, outputPath):
        """Creates a DS9-compatible .reg file with the locations of the PSF
        candidates.
        """
        serials = self.stars.keys()
        x = [self.stars[idnum]['x'] for idnum in self.stars.keys()]
        y = [self.stars[idnum]['y'] for idnum in self.stars.keys()]
        
        psfPoints = region.PointList()
        psfPoints.setFrame('image')
        psfPoints.setPoints(x, y, size=15, shapes="x", labels=serials,
                colours="red")
        psfPoints.writeTo(outputPath)
    
    # def writeWCSRegions(self, outputPath, header):
    #     """Creates a DS9 .reg file with locations of stars by their RA,Dec
    #     coordinates. The `header` is the pyfits header containing the WCS."""
    #     serials = self.stars.keys()
    #     x = [self.stars[idnum]['x'] for idnum in self.stars.keys()]
    #     y = [self.stars[idnum]['y'] for idnum in self.stars.keys()]
    #     wcs = astWCS.WCS(header, mode = "pyfits")
    #     coords = [wcs.pix2wcs(xy[0],xy[1]) for xy in zip(x,y)]
    #     ra = [coord[0] for coord in coords]
    #     dec = [coord[1] for coord in coords]
        
    #     psfPoints = region.PointList()
    #     psfPoints.setFrame('linear')
    #     psfPoints.setPoints(ra, dec, size=15, shapes="x", labels=None,
    #           colours="cyan")
    #     psfPoints.writeTo(outputPath)


class PSFModel(object):
    """Class for reading a PSF model file generated by daophot."""
    def __init__(self):
        super(PSFModel, self).__init__()
    
    def read(self, path):
        """Reads a .psf file"""
        f = open(path)
        lines = f.readlines()
        header = lines[0:2]
        f.close()
        
        line1Items = header[0].lstrip().split()
        line2Items = header[1].lstrip().split()
        
        self.psfType = line1Items[0]
        self.lutSize = int(line1Items[1])
        self.nShapeParams = int(line1Items[2])
        self.nLUT = int(line1Items[3])
        # instrumental magnitude of psf of unitary normalization
        self.mInstr = float(line1Items[5])
        # Central height, in ADU, of the analytic function which is used
        # as the first-order approximation
        self.centralHeight = float(line1Items[6])
        self.frameX0 = float(line1Items[7])
        self.frameY0 = float(line1Items[8])
        self.hwhmX = float(line2Items[0])
        self.hwhmY = float(line2Items[1])
    
    def get_prototype_mag(self):
        """Returns the magnitude of the prototype star."""
        return self.mInstr
    
    def get_seeing(self, pixelScale):
        """Returns the mean seeing (full-width at half-maximum light profile)
        of the frame in arcseconds.
        """
        meanHWHM = (self.hwhmX + self.hwhmY) / 2.
        fwhm = meanHWHM * 2.
        return fwhm * pixelScale
