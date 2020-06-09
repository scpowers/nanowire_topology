import math
import numpy as np
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
        '''
        # compute intermediate theta_s angles (see paper)
        theta_s_x = math.acos( (self.x1 - self.x_c) / self.r_w )
        theta_s_y = math.acos( (self.y1 - self.y_c) / self.r_w )
        # (experimental) compute intermediate theta_f angles
        theta_f_x = math.acos( (self.x2 - self.x_c) / self.r_w )
        theta_f_y = math.acos( (self.y2 - self.y_c) / self.r_w )

        # find theta_s from table (see paper)
        if (theta_s_x >= 0 and theta_s_x < math.pi / 2 and \
                theta_s_y >= 0 and theta_s_y < math.pi / 2):
            theta_s = theta_s_x
        elif (theta_s_x >= 0 and theta_s_x < math.pi / 2 and \
                theta_s_y >= -math.pi / 2 and theta_s_y < 0):
            theta_s = 2 * math.pi + theta_s_y
        elif (theta_s_x >= math.pi / 2 and theta_s_x < math.pi and \
                theta_s_y >= 0 and theta_s_y < math.pi / 2):
            theta_s = theta_s_x
        elif (theta_s_x >= math.pi / 2 and theta_s_x < math.pi and \
                theta_s_y >= -math.pi / 2 and theta_s_y < 0):
            theta_s = math.pi - theta_s_y
        else:
            theta_s = math.pi
            #print("entered unnatural else case on theta_s")

        # (experimental) find theta_f from table (see paper)
        if (theta_f_x >= 0 and theta_f_x < math.pi / 2 and \
                theta_f_y >= 0 and theta_f_y < math.pi / 2):
            theta_f = theta_f_x
        elif (theta_f_x >= 0 and theta_f_x < math.pi / 2 and \
                theta_f_y >= -math.pi / 2 and theta_f_y < 0):
            theta_f = 2 * math.pi + theta_f_y
        elif (theta_f_x >= math.pi / 2 and theta_f_x < math.pi and \
                theta_f_y >= 0 and theta_f_y < math.pi / 2):
            theta_f = theta_f_x
        elif (theta_f_x >= math.pi / 2 and theta_f_x < math.pi and \
                theta_f_y >= -math.pi / 2 and theta_f_y < 0):
            theta_f = math.pi - theta_f_y
        else:
            theta_f = math.pi
            #print("entered unnatural else case on theta_f")

        # switch start and finish angles, if necessary)
        if (theta_s > theta_f):
            temp = theta_s
            theta_s = theta_f
            theta_f = temp
        
        #print("Theta_s: ", theta_s)
        #print("Theta_f: ", theta_f)
        '''

        # now loop over electrodes
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

                
            
        
        



        
