import math
from scipy.ndimage import gaussian_filter
import numpy as np
import colorednoise as cn
from Electrode import Electrode
import matplotlib.pyplot as plt
from scipy.optimize import minimize 
import copy
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
    def setScaledPoints(self, length, electrodes):
        # make length an attribute for later optimization
        self.length = length
        # define coords of the point about which the polynomial
        # is to be rotated (for later)
        self.x_r = length / 2
        self.y_r = self.x_r

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

        # do secondary transformation of y values
        tempMin = min(self.scaledPoints)
        tempMax = max(self.scaledPoints)
        for i in range(self.numPoints):
            v = np.random.uniform(-1 * tempMin, length - tempMax)
            self.scaledPoints[i] = self.scaledPoints[i] + v

        # feed y-values through Gaussian filter function (per paper)
        self.scaledPoints = gaussian_filter(self.scaledPoints, sigma=3)

        # fit quintic polynomial to data
        # make it a self attribute for reference in constraint function
        self.poly_coeff = np.polyfit(self.scaledXPoints, self.scaledPoints, 5)

        # store the function output at each x in place of the scaled y-values
        for i in range(self.numPoints):
            self.scaledPoints[i] = self.quintic(self.scaledXPoints[i])
            
        # add shift to curve to create better distribution
        self.y_rand = np.random.uniform(-1 * length / 2, 0)
        self.scaledPoints = self.scaledPoints  + self.y_rand

            
        # rotate the fitted function about an arbitrary point by a random theta
        # they chose the rotation point (l/2, l/2)
        self.theta = np.random.uniform(0, 2 * math.pi)
        x_r = length / 2
        y_r = length / 2
        self.x_rotated, self.y_rotated = self.rotateXY(self.scaledXPoints, 
                self.scaledPoints, x_r, y_r, self.theta)

        # check for electrode connectivity while you have access to the
        # coefficient array from the polynomial fit
        self.checkIfConnected(electrodes)
                
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
    def quintic(self, x_in):
        return (self.poly_coeff[0] * (x_in ** 5) + self.poly_coeff[1] * \
                (x_in ** 4) + self.poly_coeff[2] * \
                (x_in ** 3) + self.poly_coeff[3] * (x_in ** 2) + \
                self.poly_coeff[4] * (x_in) + self.poly_coeff[5])

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
        
    # override inherited method checking for connected electrodes
    def checkIfConnected(self, electrodes):
        # loop over the list of electrodes
        for electrode in electrodes:
            # get coordinates of the electrode center
            [self.x_e, self.y_e] = electrode.getCenter()
            # get electrode radius
            rad_e = electrode.getRadius()

            # get roots of polynomial
            rootSol = minimize(self.d_squared, method='Nelder-Mead', x0=0) 
            roots = rootSol.x
            if len(roots) == 0:
                print("no real roots found!")
            
            # roots contain x values...we need to run them
            # back through the rotation process and then
            # we can get the distance
            # check output
            for root in roots:
                # get unrotated y
                y_temp = self.quintic(root)
                y_temp = y_temp + self.y_rand

                # get rotated x and y
                x_rot, y_rot = self.rotateXY([root], [y_temp],
                        self.x_r, self.y_r, self.theta)

                # compute distance between electrode center
                # and rotated x and y from root
                dist = math.sqrt( (x_rot[0] - self.x_e) ** 2 + \
                        (y_rot[0] - self.y_e) ** 2 )

                # check condition
                if dist <= rad_e:
                    self.connectedElectrodes.append(electrode)
                    # if this root intersects this electrode,
                    # no need to check other roots for this 
                    # same electrode
                    break
                    
    # now the function to minimize (square of distance)
    def d_squared(self, x):
        # get unrotated y coordinate
        y_temp = self.quintic(x)
        y_temp = y_temp + self.y_rand
        # get rotated points
        rx_temp, ry_temp = self.rotateXY([x], [y_temp], 
                self.x_r, self.y_r, self.theta)
        return (rx_temp[0] - self.x_e) ** 2 + (ry_temp[0] - self.y_e) ** 2
        
    
    # method for printing useful info about the wire
    def __str__(self):
        temp = []
        for electrode in self.connectedElectrodes:
            temp.append(electrode.getIndex())
        
        ret = "Wire details:\n\tEndpoint 1: (" + str(self.x_rotated[0]) + ", " + \
                str(self.y_rotated[0]) + ")" \
            + "\n\tEndpoint 2: (" + str(self.x_rotated[self.numPoints-1]) + ", " + \
            str(self.y_rotated[self.numPoints-1]) + ")" + \
            "\n\tConnected Electrodes: " + str(temp) + "\n\tIndex: " + str(self.index)
        return ret

        
