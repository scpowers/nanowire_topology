import numpy as np
import matplotlib.pyplot as plt
from GeometricModel import GeometricModel
from Wire import StraightWire

# StraightWireModel class: defines a type of GeometricModel that uses straight
# lines to model the memristive nanowires
class StraightWireModel(GeometricModel):

    # inherit the __init__() method from parent class, no added parameters
    def __init__(self, density, r_e, alpha, num_e):
        GeometricModel.__init__(self, density, r_e, alpha, num_e)
        
        # initialize wire list
        self.Wires = []
    
    # generate wire models
    def generateWires(self):
        for wire in range(self.num_w):
            
            # draw side values (can't be the same for straight wires)
            q = np.random.choice( range(4), 2, replace=False )

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

            # create straight wire object
            self.Wires.append(StraightWire(coords)) 
            # update list of connected electrodes
            self.Wires[wire].checkIfConnected(self.ElectrodeList)
            print(self.Wires[-1])

    def plotModel(self):

        # initialize figure
        plt.figure(1)
        ax = plt.gca()
        
        # loop over electrodes
        for electrode in range(self.num_e):
            coords = self.ElectrodeList[electrode].getCenter()
            rad = self.ElectrodeList[electrode].getRadius()
            circle = plt.Circle((coords[0], coords[1]), rad, color='r')
            ax.add_artist(circle)

        for wire in range(self.num_w):
            coords = self.Wires[wire].getEndpoints()
            x1 = coords[0][0]
            y1 = coords[0][1]
            x2 = coords[1][0]
            y2 = coords[1][1]
            if ( (x1 > 0 and x1 < self.length and y1 > 0 and y1 < self.length) or \
                    (x2 > 0 and x2 < self.length and y2 > 0 and y2 < self.length) ):
                print("found a problem: " + str(x1) + ", " + str(y1) + \
                    "\n" + str(x2) + ", " + str(y2))
            plt.plot([x1, x2], [y1, y2], color='blue')

        plt.savefig('plot.png')
        plt.xlim(0, self.length)
        plt.ylim(0, self.length)
        plt.show()

        
        

