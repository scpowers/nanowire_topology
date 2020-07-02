import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from GeometricModel import GeometricModel
from StraightWireModel import StraightWireModel
from Wire import *
from helpers import generateBipartiteGraph

# ArcWireModel class: defines a model that uses arcs of random curvature
# and centerpoint to define wires (more realistic than straight lines)
class ArcWireModel(StraightWireModel):

    
    # inherit the __init__() method from parent class, no added parameters
    def __init__(self, density, r_e, alpha, num_e):
        GeometricModel.__init__(self, density, r_e, alpha, num_e)
        
        # generate wires
        self.generateWires()
    
    # generate wire models
    def generateWires(self):
        for wire in range(self.num_w):
            
            # draw side values (can be the same for arced wires)
            q = np.random.randint( 4, size=2 )

            # draw coordinates from uniform distribution
            coords = np.random.uniform(0, self.length, (2,2))
            
            # update coordinates based on table in Pantone et al
            for point in range(2):
                if (q[point] == 0):
                    coords[point, 1] = 0
                elif (q[point] == 1):
                    coords[point, 0] = 0
                elif (q[point] == 2):
                    coords[point, 1] = self.length
                else:
                    coords[point, 0] = self.length

            # create arc wire object
            self.Wires.append(ArcWire(coords, wire)) 
            # update list of connected electrodes
            self.Wires[wire].checkIfConnected(self.ElectrodeList)
            #print(self.Wires[-1])
            
        
    # override plot geometric model (electrodes and wires)
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
            [x_c, y_c, r] = self.Wires[wire].getPlottingInfo()
            circle = plt.Circle((x_c, y_c), r, color='blue', fill=False)
            ax.add_artist(circle)
            #print("x_c, y_c, r: ", x_c, " ", y_c, " ", r)

        plt.xlim(0, self.length)
        plt.ylim(0, self.length)
        plt.savefig('GeometricModelPlot.png')
        plt.close()
        #plt.show()

