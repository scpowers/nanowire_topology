import math
from scipy.ndimage import gaussian_filter
import numpy as np
import colorednoise as cn
from Electrode import Electrode
import matplotlib.pyplot as plt
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
        # scale them to the correct range
        for i in range(self.numPoints):
            self.scaledPoints[i] = (length * (self.scaledPoints[i] - tempMin)) / (tempMax - tempMin)

        # set scaled x values while you're at it
        # note that I had to extend the x range to ensure that the ends of the wire
        # did not end up in the work area after rotation
        bufferLength = (length * (math.sqrt(2) - 1)) / 2
        self.scaledXPoints = np.linspace(-1 * bufferLength, length + bufferLength, 
                len(self.scaledPoints))

        #print("after initial scaling: y min/max = ", min(self.scaledPoints), " ", 
        #        max(self.scaledPoints))
        """
        plt.plot(self.scaledXPoints, self.scaledPoints)
        plt.xlim(-1 * bufferLength, length + bufferLength)
        plt.ylim(-1 * bufferLength, length + bufferLength)
        plt.show()
        """
        
        # do secondary transformation of y values
        for i in range(self.numPoints):
            r = np.random.uniform(0, 1)
            self.scaledPoints[i] = self.scaledPoints[i] * (r ** 1) # from paper

        #print("after multiplying by r^3: y min/max = ", min(self.scaledPoints), " ", 
        #        max(self.scaledPoints))
        """
        plt.plot(self.scaledXPoints, self.scaledPoints)
        plt.xlim(0, length)
        plt.ylim(0, length)
        plt.show()
        """
        
        # do tertiary transformation of y values
        tempMin = min(self.scaledPoints)
        tempMax = max(self.scaledPoints)
        # temp copy for comparison later
        tempCopy = self.scaledPoints.copy()
        for i in range(self.numPoints):
            #print("point before: ", point)
            v = np.random.uniform(-1 * tempMin, length - tempMax)
            #print("v = ", v)
            self.scaledPoints[i] = self.scaledPoints[i] + v
            #print("point after: ", point)
        #print("did anything change? diff vector: ", tempCopy - self.scaledPoints)

        #print("after adding v: y min/max = ", min(self.scaledPoints), " ", 
        #        max(self.scaledPoints))
        
        """
        plt.plot(self.scaledXPoints, self.scaledPoints)
        plt.xlim(0, length)
        plt.ylim(0, length)
        plt.show()
        """
        
        # feed y-values through Gaussian filter function (per paper)
        # guess 1 for SD..not listed in paper
        self.scaledPoints = gaussian_filter(self.scaledPoints, sigma=3)

        # interested to see what the range of y values is after filtering
        #print("after filtering: y min/max = ", min(self.scaledPoints), " ", 
        #        max(self.scaledPoints))

        # fit quintic polynomial to data
        poly_coeff = np.polyfit(self.scaledXPoints, self.scaledPoints, 5)

        """
        plt.plot(self.scaledXPoints, self.scaledPoints)
        plt.xlim(0, length)
        plt.ylim(0, length)
        plt.show()
        """
        
        # store the function output at each x in place of the scaled y-values
        for i in range(self.numPoints):
            self.scaledPoints[i] = self.quintic(poly_coeff, self.scaledXPoints[i])
            
        # rotate the fitted function about an arbitrary point by a random theta
        # they chose the rotation point (l/2, l/2)
        theta = np.random.uniform(0, 2 * math.pi)
        x_r = length / 2
        y_r = length / 2
        #self.x_rotated, self.y_rotated = self.rotateXY(self.scaledXPoints, 
        #        self.scaledPoints, x_r, y_r, theta)
        # don't rotate for now
        self.x_rotated = self.scaledXPoints
        self.y_rotated = self.scaledPoints
        #print("After rotation:")
        #print("x min/max: ", min(self.x_rotated), " ", max(self.x_rotated))
        #print("y min/max: ", min(self.y_rotated), " ", max(self.y_rotated))

        """
        plt.plot(self.x_rotated, self.y_rotated)
        plt.xlim(0, length)
        plt.ylim(0, length)
        plt.show()
        """
        


        # random test - remove later
        """
        t_test = np.arange(0, math.pi, math.pi/30)
        x_test = np.zeros(len(t_test))
        y_test = np.zeros(len(t_test))
        for t in range(len(t_test)):
            x_test[t] = 1.5 + 1.5 * math.cos(t_test[t])
            y_test[t] = 1.5 * math.sin(t_test[t])
        print("min/max x: ", min(x_test), " ", max(x_test))
        print("min/max y: ", min(y_test), " ", max(y_test))
        plt.plot(x_test, y_test)
        plt.xlim(0, 3)
        plt.ylim(0, 3)
        plt.show()

        # rotate the random test
        x_rot, y_rot = self.rotateXY(x_test, y_test, x_r, y_r, 0.01)
        plt.plot(x_rot, y_rot)
        plt.xlim(0, 3)
        plt.ylim(0, 3)
        plt.show()
        """ 

    # helper function for rotating points about another point by an angle
    def rotateXY(self, x, y, x_r, y_r, theta):
        # initialize arrays of the same size (will get overwritten)
        r_x = x.copy()
        r_y = y.copy()
        for i in range(len(x)):
            # apply rotation matrix operation to x coordinate
            r_x[i] = (x[i] - x_r) * math.cos(theta) - \
                    (y[i] - y_r) * math.sin(theta) + x_r
            
            # apply rotation matrix operation to y coordinate
            r_y[i] = (x[i] - x_r) * math.sin(theta) + \
                    (y[i] - y_r) * math.cos(theta) + y_r
        return r_x, r_y
        
    # helper function for the quintic polynomial
    def quintic(self, coeff_vec, x_in):
        return (coeff_vec[0] * (x_in ** 5) + coeff_vec[1] * (x_in ** 4) + coeff_vec[2] * \
                (x_in ** 3) + coeff_vec[3] * (x_in ** 2) + coeff_vec[4] * (x_in) + coeff_vec[5])

    # helper function for the quintic polynomial
    def sept(self, coeff_vec, x_in):
        return (coeff_vec[0] * (x_in ** 7) + coeff_vec[1] * (x_in ** 6) + coeff_vec[2] * \
                (x_in ** 5) + coeff_vec[3] * (x_in ** 4) + coeff_vec[4] * (x_in ** 3) + \
                coeff_vec[5] * (x_in ** 2) + coeff_vec[6] * (x_in) + coeff_vec[7])
        
    # getter for final rotated x points
    def getRotatedXPoints(self):
        return self.x_rotated

    # getter for final rotated y points
    def getRotatedYPoints(self):
        return self.y_rotated
        
