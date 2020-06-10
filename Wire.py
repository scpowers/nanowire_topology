import math
from scipy.ndimage import gaussian_filter
import numpy as np
import colorednoise as cn
from Electrode import Electrode
# Wire class - defines attributes of various wire models

class StraightWire:

    def __init__(self, coords, index):
        # read in endpoint coordinates
        self.x1 = coords[0,0] # x coordinate of first endpoint
        self.y1 = coords[0,1] # y coordinate of first endpoint
        self.x2 = coords[1,0] # x coordinate of second endpoint
        self.y2 = coords[1,1] # y coordinate of second endpoint
        self.index = index
        
        # initialize array of electrodes connected by this wire
        self.connectedElectrodes = []

    # getter method for endpoints
    def getEndpoints(self):
        # row 1 is endpoint 1, row 2 is endpoint 2
        return [[self.x1, self.y1],
                [self.x2, self.y2]]

    # getter method for connected electrodes
    def getConnectedElectrodes(self):
        return self.connectedElectrodes

    # getter method for index
    def getIndex(self):
        return self.index
        
    # method for printing useful info about the wire
    def __str__(self):
        temp = []
        for electrode in self.connectedElectrodes:
            temp.append(electrode.getIndex())
        
        ret = "Wire details:\n\tEndpoint 1: (" + str(self.x1) + ", " + str(self.y1) + ")" \
            + "\n\tEndpoint 2: (" + str(self.x2) + ", " + str(self.y2) + ")" + \
            "\n\tConnected Electrodes: " + str(temp) + "\n\tIndex: " + str(self.index)
        return ret

    # check if electrodes in a list of electrodes touch this wire
    def checkIfConnected(self, electrodes):
        for electrode in electrodes:
            
            # for comparison purposes
            rad = electrode.getRadius()

            # get center coordinates of the electrode
            c = electrode.getCenter()

            # see paper for equation
            num = abs( c[0] * (self.y2 - self.y1) - c[1] * (self.x2 - self.x1) + \
                    self.x2 * self.y1 - self.y2 * self.x1 )
            denom = math.sqrt( (self.x2 - self.x1) ** 2 + (self.y2 - self.y1) \
                    ** 2 )
            D_s = num / denom

            # if it's overlapping, add this electrode to the list
            if (D_s <= rad):
                self.connectedElectrodes.append(electrode)


# class describing the arced wire model
class ArcWire(StraightWire):

    # call inherited constructor
    def __init__(self, coords, index):
        # call inherited constructor
        StraightWire.__init__(self, coords, index)

        # derived attributes
        self.chooseRadius()
        self.findCenter()


    # define helper method for randomly choosing wire arc radius
    def chooseRadius(self):
        Euclid_dist = math.sqrt( (self.x2 - self.x1) ** 2 + (self.y2 - \
                self.y1) ** 2 )
        lowerBound = Euclid_dist / 2 # from paper
        upperBound = 3 * Euclid_dist # from paper
        self.r_w = np.random.uniform(lowerBound, upperBound)

    # define helper method for computing center coordinates
    def findCenter(self):

        # intermediate variables
        x_a = (self.x2 - self.x1) / 2
        y_a = (self.y2 - self.y1) / 2
        a = math.sqrt( x_a ** 2 + y_a ** 2 )
        b = math.sqrt( self.r_w ** 2 - a ** 2 )

        # randomly choose one of the two possible centers
        x_rand = np.random.randint(2) # binary choice of 0 or 1
        #print("random value: ", x_rand)

        # assign values based on random choice above
        if (x_rand == 0):
            self.x_c = self.x1 + x_a + (b * y_a / a)
            self.y_c = self.y1 + y_a - (b * x_a / a)
        else:
            self.x_c = self.x1 + x_a - (b * y_a / a)
            self.y_c = self.y1 + y_a + (b * x_a / a)

    # return center and radius for plotting
    def getPlottingInfo(self):
        return [self.x_c, self.y_c, self.r_w]

    # override inherited method checking for connected electrodes
    def checkIfConnected(self, electrodes):

        # loop over electrodes
        for electrode in electrodes: 
            # for comparison purposes
            rad = electrode.getRadius()

            # get center coordinates of the electrode
            [x_e, y_e] = electrode.getCenter()

            # compute distance between arc center and electrode center
            dist = math.sqrt( (x_e - self.x_c) ** 2 + (y_e - self.y_c) ** 2 )
            
            # evaluate condition
            if ( dist >= self.r_w - rad and dist <= self.r_w + rad ):
                self.connectedElectrodes.append(electrode)

# class describing the pink noise wire model
class PinkNoiseWire(StraightWire):

    # define new constructor
    def __init__(self, index, numPoints):
            
        # set wire index
        self.index = index
        
        # initialize array of electrodes connected by this wire
        self.connectedElectrodes = []
        
        # new attributes
        self.numPoints = numPoints

        # called methods for derived attributes
        # add them here...these calls must define
        # all behavior of the wire
        self.generatePinkNoise()

    # method for generating pink noise points for this wire
    def generatePinkNoise(self):
        # using colorednoise Python package to generate pink noise
        # beta = 1 for pink noise
        self.rawPoints = cn.powerlaw_psd_gaussian(1, self.numPoints)

    # setter method for setting rescale points
    def setScaledPoints(self, length):

        # rescale raw noise (y-values)
        tempMin = min(self.rawPoints)
        tempMax = max(self.rawPoints)
        # copy the raw points
        self.scaledPoints = self.rawPoints[:]
        for point in self.scaledPoints:
            point = (length * (point - tempMin)) / (tempMax - tempMin)

        # set scaled x values while you're at it
        self.scaledXPoints = np.linspace(0, length, len(self.scaledPoints))

        # do secondary transformation of y values
        for point in self.scaledPoints:
            r = np.random.uniform()
            point = point * (r ** 3) # from paper

        # do tertiary transformation of y values
        tempMin = min(self.scaledPoints)
        tempMax = max(self.scaledPoints)
        for point in self.scaledPoints:
            v = np.random.uniform(tempMin, tempMax)
            point = point + v

        # feed y-values through Gaussian filter function (per paper)
        # guess 1 for SD..not listed in paper
        self.scaledPoints = gaussian_filter(self.scaledPoints, sigma=1)

        # fit quintic polynomial to data
        poly_coeff = np.polyfit(self.scaledXPoints, self.scaledPoints, 5)

        # store the function output at each x in place of the scaled y-values
        for i in range(self.numPoints):
            self.scaledPoints[i] = self.quintic(poly_coeff, self.scaledXPoints[i])
            
        # rotate the fitted function about an arbitrary point by a random theta
        # they chose the rotation point (l/2, l/2)
        theta = np.random.uniform(0, 2 * math.pi)
        x_r = length / 2
        y_r = length / 2
        self.x_rotated = np.zeros(self.numPoints)
        self.y_rotated = np.zeros(self.numPoints)
        for i in range(self.numPoints):
            self.x_rotated[i] = (self.scaledXPoints[i] - x_r) * math.cos(theta) - \
                    (self.scaledPoints[i] - y_r) * math.sin(theta) + x_r
            self.y_rotated[i] = (self.scaledXPoints[i] - x_r) * math.sin(theta) + \
                    (self.scaledPoints[i] - y_r) * math.cos(theta) + y_r

        
    # helper function for the quintic polynomial
    def quintic(self, coeff_vec, x_in):
        return coeff_vec[0] * (x_in ** 5) + coeff_vec[1] * (x_in ** 4) + coeff_vec[2] * \
                (x_in ** 3) + coeff_vec[3] * (x_in ** 2) + coeff_vec[4] * (x_in) + coeff_vec[5]

    # getter for final rotated x points
    def getRotatedXPoints(self):
        return self.x_rotated

    # getter for final rotated y points
    def getRotatedYPoints(self):
        return self.y_rotated
        
