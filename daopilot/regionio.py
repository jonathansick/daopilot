"""region module. Makes DS9-compatible region files and deals with converting
regions/points in world coordinates to image coordiantes, and vice versa.

See DS9 formatting information at:
http://hea-www.harvard.edu/RD/ds9/ref/region.html#RegionDescriptions
"""

import os
import re


class PointList(object):
    """docstring for PointList"""
    def __init__(self):
        super(PointList, self).__init__()
        self.frame = 'fk5'
        self.x = None
        self.y = None
        self.labels = None
        self.colours = None
        self.shapes = "circle"  # or box, diamond, etc,
        self.size = 4

    def set_frame(self, frameType):
        """The coordinate frame can either be 'image' or 'fk5'."""
        self.frame = frameType

    def set_points(self, x, y, shapes="circle", labels=None, colours=None,
            size=4):
        """Sets point data (x,y) or (ra,dec) and optional labels as a list
        of python strings.
        
        white
        black
        red
        green
        blue
        cyan
        magenta
        yellow
        """
        self.x = x
        self.y = y
        self.size = size
        n = len(self.x)
        
        if type(shapes) is list or type(shapes) is tuple:
            self.shapes = shapes
        else:
            self.shapes = [shapes] * n
        
        if labels is None or type(labels) is list or type(labels) is tuple:
            self.labels = labels
        else:
            self.labels = [labels] * n
        
        if colours is None or type(colours) is list or type(colours) is tuple:
            self.colours = colours
        else:
            self.colours = [colours] * n

    def make_lines(self):
        """Returns text of the point data in DS9 region format."""
        n = len(self.x)
        lineList = []
        # print self.shapes
        for i in xrange(n):
            line = "point(%f,%f) # point=%s %i" % (self.x[i], self.y[i],
                self.shapes[i], self.size)
            
            if self.labels is not None:
                line = " ".join((line, "text = {%s}" % self.labels[i]))
            if self.colours is not None:
                line = " ".join((line, "color = %s" % self.colours[i]))
            lineList.append(line)
        
        return "\n".join(lineList)
    
    def write_to(self, outputPath):
        """Writes the region data to the *outputPath*."""
        text = "\n".join((self.frame, self.makeLines()))
        if os.path.exists(os.path.dirname(outputPath)) is False:
            os.makedirs(os.path.dirname(outputPath))
        if os.path.exists(outputPath):
            os.remove(outputPath)
        f = open(outputPath, 'w')
        f.write(text)
        f.close()


class BoxList(object):
    """Reads and writes lists of boxes in DS9 .reg format."""
    def __init__(self):
        super(BoxList, self).__init__()
        self.regions = None
        self.frame = 'fk5'
    
    def read(self, path):
        """At this point, we need the region to have fk5 coordinates in
        hexagesimal format.
        """
        self.regions = []
        f = open(path, 'rU')
        for line in f:
            if line.startswith("box") is False:
                continue
            
            # parse the coordinates
            matches = re.match(r'box\(([0-9][0-9]):([0-9][0-9]):([0-9\.]*),([\+\-0-9]*):([0-9]*):([0-9\.]*)', line)
            if matches is None:
                continue
            
            items = matches.groups()
            raH = int(items[0])
            raM = int(items[1])
            raS = float(items[2])
            ra = (raH + raM / 60. + raS / 3600.) * 15.
            decD = int(items[3])
            decM = int(items[4])
            decS = float(items[5])
            dec = decD + decM / 60. + decS / 3600.
            
            region = {'raH': raH, 'raM': raM, 'raS': raS, 'ra': ra,
                'decD': decD, 'decM': decM, 'decS': decS, 'dec': dec}
            
            # parse the text label
            matches = re.match(r'.*text={([0-9a-zA-Z]*)}', line)
            if matches:
                text = matches.groups()[0]
                region['text'] = text
            
            self.regions.append(region)
        f.close()
    
    def get_points(self):
        """Returns the list of box dictionaries"""
        return self.regions
    
    def set_points(self, x, y, xSize, ySize, text, frame='fk5'):
        """docstring for setPoints"""
        self.frame = frame
        self.regions = []
        for i in xrange(len(x)):
            region = {'x': x[i], 'y': y[i], 'x_size': xSize[i],
                    'y_size': ySize[i], 'text': text[i]}
            self.regions.append(region)

    def write_to(self, path):
        """docstring for writeTo"""
        lineList = [self.frame]
        for reg in self.regions:
            line = "box(%s,%s,%s,%s) # text={%s}" % \
                (reg['x'], reg['y'], reg['x_size'], reg['y_size'], reg['text'])
            lineList.append(line)
        if os.path.exists(path):
            os.remove(path)
        f = open(path, 'w')
        f.write("\n".join(lineList))
        f.close()
