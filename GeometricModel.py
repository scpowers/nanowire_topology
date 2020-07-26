import math
import numpy as np
from Electrode import Electrode

# class GeometricModel: builds a geometric model of the memristive nanowire array
class GeometricModel:

    def __init__(self, density, r_e, alpha, num_e, spacing=False):
        # user-specified model parameters
        self.density = density # wire density constant
        self.r_e = r_e # electrode radius
        self.alpha = alpha # distance between electrodes on square
        self.num_e = num_e # number of electrodes (perfect square)
        self.spacing = spacing # boolean for clustering electrodes

        # derived parameters
        self.numRows = int(math.sqrt(self.num_e)) # number of rows and columns of electrodes
        self.num_w = self.density * self.numRows # number of wires
        self.length = self.alpha * (1 + self.numRows) # grid side length
        self.ElectrodeList = [] # list of Electrode objects
        self.Wires = []

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
                
                # determine the index
                index = col + (row * self.numRows)
                
                # if you want to separate clusters of electrodes, apply modifications
                if (self.spacing == True):
                    x_temp = index % math.sqrt(self.num_e) 
                    y_temp = index // math.sqrt(self.num_e)
                    numDivs = int(math.log2(math.sqrt(self.num_e) // 2))

                    for i in range(1,numDivs+1):
                        tmp = 2 ** i
                        
                        x_buffer = (x_temp // tmp) * (tmp * self.alpha)

                        y_buffer = (y_temp // tmp) * (tmp * self.alpha)
                        
                        x += x_buffer
                        y += y_buffer
                    
                # create new Electrode object at the right coordinates (constructor) 
                tempElectrode = Electrode(self.r_e, x, y, index)
                self.ElectrodeList.append(tempElectrode)
                #print(self.ElectrodeList[-1])
                 
            
        
