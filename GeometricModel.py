import math
import numpy as np
from Electrode import Electrode

# class GeometricModel: builds a geometric model of the memristive nanowire array
class GeometricModel:

    def __init__(self, density, r_e, alpha, num_e):
        # user-specified model parameters
        self.density = density # wire density constant
        self.r_e = r_e # electrode radius
        self.alpha = alpha # distance between electrodes on square
        self.num_e = num_e # number of electrodes (perfect square)

        # derived parameters
        self.numRows = int(math.sqrt(self.num_e)) # number of rows and columns of electrodes
        self.num_w = self.density * self.numRows # number of wires
        self.length = self.alpha * (1 + self.numRows) # grid side length
        self.ElectrodeList = [] # list of Electrode objects

        # generate electrodes
        self.initializeElectrodes()

    # method for initializing the electrodes (universal across all sub-models)
    def initializeElectrodes(self):
        edgeBuffer = ( self.length - ( self.numRows - 1) * self.alpha ) / 2

        for row in range(self.numRows):
            for col in range(self.numRows):
                # determine x and y for this electrode
                x = edgeBuffer + ( col * self.alpha )
                y = edgeBuffer + ( row * self.alpha )
                # create new Electrode object at the right coordinates (constructor) 
                index = col + (row * self.numRows)
                tempElectrode = Electrode(self.r_e, x, y, index)
                self.ElectrodeList.append(tempElectrode)
                #print(self.ElectrodeList[-1])
                 
        
