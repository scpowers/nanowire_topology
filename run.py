#!/usr/bin/env python3
import sys
import os
sys.path.append(os.getcwd())

# run.py - file to actually run
from GeometricModel import GeometricModel
from StraightWireModel import StraightWireModel
from Electrode import Electrode

# model parameters
r_e = 0.4 # electrode radius
alpha = 1 # distance between electrode centers
density = 30 # wire density constant
num_e = 4 # number of electrodes (perfect square)

model = StraightWireModel(density, r_e, alpha, num_e)
model.initializeElectrodes()
model.generateWires()
model.plotModel()


