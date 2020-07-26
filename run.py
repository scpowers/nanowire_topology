#!/usr/bin/env python3
import sys
import os
sys.path.append(os.getcwd())

# run.py - file to actually run
from GeometricModel import GeometricModel
from StraightWireModel import StraightWireModel
from PinkNoiseModel import PinkNoiseModel
from Electrode import Electrode
from helpers import *
import networkx as nx
from random import sample
import matplotlib.pyplot as plt

# model parameters
r_e = 0.4 # electrode radius
alpha = 1 # distance between electrode centers
density = 30 # wire density constant

# specific for pink noise model
numPointsPerWire = 201 # from paper

# sweep parameters for the square root of the number of electrodes
num_e_min = 3
num_e_max = 40

# initialize array for small world coefficient parameters
sigma_array = []

#"""
# begin sweep of square root of electrode count:
for root_num_e in range(num_e_min, num_e_max + 1):
    
    print("working on root_num_e = ", root_num_e)
    # initialize geometric model
    model = PinkNoiseModel(density, r_e, alpha, root_num_e ** 2, numPointsPerWire, spacing=False)
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
#"""

"""
# create large network first
model = PinkNoiseModel(density, r_e, alpha, 128 ** 2, numPointsPerWire, spacing=True)
# plot model
model.plotModel()
# generate equivalent bipartite graph
equivalentGraph = model.generateGraph()
#nx.draw(equivalentGraph)
#plt.show()

# sweep chemical distance and store N(r) in a vector
r_vec = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
N_vec = []


random_nodes = sample(list(equivalentGraph.nodes()), len(r_vec))

for i in range(len(r_vec)):
    r = r_vec[i]
    seed = random_nodes[i]
    ego_net = nx.ego_graph(equivalentGraph, seed, radius=r)
    N_vec.append(len(ego_net.nodes()))
    print("for r = ", r_vec[i], ": N(r) = ", N_vec[i])

print(r_vec)
print(N_vec)
plt.loglog(r_vec, N_vec)
plt.xlabel('r')
plt.ylabel('N(r)')
plt.savefig('N_vs_r_plot.png')
#"""
