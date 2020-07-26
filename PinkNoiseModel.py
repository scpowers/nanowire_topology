import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from GeometricModel import GeometricModel
from Wire import *
from helpers import generateBipartiteGraph
from StraightWireModel import StraightWireModel

# PinkNoiseModel class: defines a model that uses points generated
# via pink noise to define random 5th-order polynomials that are
# rotated about a random point to produce realistic nanowires
class PinkNoiseModel(StraightWireModel):

    # call inherited constructor
    def __init__(self, density, r_e, alpha, num_e, numPointsPerWire, spacing=False):
        GeometricModel.__init__(self, density, r_e, alpha, num_e, spacing)

        # new attribute: number of points per wire
        self.numPoints = numPointsPerWire

        # rescale limits of boundary box given spacing boolean
        self.lim = self.length
        if (self.spacing == True):
            limCopy = self.lim
            numDivs = int(math.log2(math.sqrt(self.num_e) // 2))

            for i in range(1,numDivs+1):
                tmp = 2 ** i
                self.lim += ((math.sqrt(self.num_e) // tmp ) - 1) * (tmp * self.alpha)

        
        #if (self.spacing == True):
        #    self.lim += ((self.lim // 4) - 1) * (2 * self.alpha) + \
        #            ((self.lim // 8) - 1) * (4 * self.alpha) + \
        #            ((self.lim // 16) - 1) * (8 * self.alpha)
                    
        # generate wires
        self.generateWires()


    # generate wire models
    def generateWires(self):

        for wire in range(self.num_w):

            # initialize new wire
            newWire = PinkNoiseWire(wire, self.numPoints) 

            # set scaled points (both x and y) and apply another transformation
            # on the y-values (see paper)
            newWire.setScaledPoints(self.lim, self.ElectrodeList)

            # add wire to list
            self.Wires.append(newWire)

            #print(self.Wires[wire])

    # override plot method for geometric model (electrodes and wires)
    def plotModel(self):

        # initialize figure
        plt.figure(1)
        ax = plt.gca()
        
        # loop over electrodes
        for electrode in range(self.num_e):
            coords = self.ElectrodeList[electrode].getCenter()
            rad = self.ElectrodeList[electrode].getRadius()
            # add circle object
            circle = plt.Circle((coords[0], coords[1]), rad, color='r')
            ax.add_artist(circle)

        # loop over wires
        for wire in range(self.num_w):
            x = self.Wires[wire].getRotatedXPoints()
            y = self.Wires[wire].getRotatedYPoints()
            plt.plot(x, y, color='blue')

        # generate plot
        
        plt.xlim(0, self.lim)
        plt.ylim(0, self.lim)
        plt.savefig('GeometricModelPlot.png')
        #plt.close()
        #plt.show()
        
        


        
        

