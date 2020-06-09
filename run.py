#!/usr/bin/env python3
import sys
import os
sys.path.append(os.getcwd())

# run.py - file to actually run
from GeometricModel import GeometricModel
from StraightWireModel import *
from Electrode import Electrode
from helpers import *

# model parameters
r_e = 0.4 # electrode radius
alpha = 1 # distance between electrode centers
density = 30 # wire density constant

# sweep parameters for the square root of the number of electrodes
num_e_min = 3
num_e_max = 15

# initialize array for small world coefficient parameters
sigma_array = []

# begin sweep of square root of electrode count:
for root_num_e in range(num_e_min, num_e_max + 1):
    
    print("working on root_num_e = ", root_num_e)
    # initialize geometric model
    model = ArcWireModel(density, r_e, alpha, root_num_e ** 2)
    # plot geometric model
    model.plotModel()
    # generate equivalent bipartite graph
    [equivalentGraph, e_nodes] = model.generateEquivalentBipartiteGraph()
    # generate random bipartite graph with the same node and edge counts
    [randomGraph, rand_e_nodes] = generateRandomBipartiteGraph(model.num_e, 
            model.numValidWires, model.numValidEdges)
    # get average shortest path length and square clustering coefficient for both graphs
    [L, C] = analyzeBipartiteGraph(equivalentGraph, e_nodes)
    [L_r, C_r] = analyzeBipartiteGraph(randomGraph, rand_e_nodes)

    # compute small world coefficient
    sigma_array.append(computeSmallWorldCoefficient(L, C, L_r, C_r))

# plot small world coefficient versus square root of electrode count
plotSmallWorldSweep(sigma_array, range(num_e_min, num_e_max + 1))


