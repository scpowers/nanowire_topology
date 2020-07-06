import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from GeometricModel import GeometricModel
from Wire import *
from helpers import generateBipartiteGraph

# StraightWireModel class: defines a type of GeometricModel that uses straight
# lines to model the memristive nanowires
class StraightWireModel(GeometricModel):

    # inherit the __init__() method from parent class, no added parameters
    def __init__(self, density, r_e, alpha, num_e):
        GeometricModel.__init__(self, density, r_e, alpha, num_e)
        
        # generate wires
        self.generateWires()
    
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
            self.Wires.append(StraightWire(coords, wire)) 
            # update list of connected electrodes
            self.Wires[wire].checkIfConnected(self.ElectrodeList)
            #print(self.Wires[-1])

    # plot geometric model (electrodes and wires)
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
            coords = self.Wires[wire].getEndpoints()
            x1 = coords[0][0]
            y1 = coords[0][1]
            x2 = coords[1][0]
            y2 = coords[1][1]
            # plot lines based on two coordinates
            plt.plot([x1, x2], [y1, y2], color='blue')

        plt.xlim(0, self.length)
        plt.ylim(0, self.length)
        plt.savefig('GeometricModelPlot.png')
        #plt.close()
        plt.show()

    # generate equivalent bipartite graph
    def generateEquivalentBipartiteGraph(self):

        # first, generate string labels for electrode nodes (need to be distinct from wire nodes)
        e_nodes = []
        for electrode in self.ElectrodeList:
            e_nodes.append("e" + str(electrode.getIndex()))
        #print(e_nodes)

        # next, generate string labels for wire nodes
        # and also generate edges simultaneously
        w_nodes = [] # list of wire node labels
        edges_list = [] # list of edges
        for wire in self.Wires:
            # only care if the wire connects 2 or more electrodes
            numConnected = len(wire.getConnectedElectrodes())
            if numConnected >= 2:
                w_nodes.append("w" + str(wire.getIndex()))

                # build edges
                tempList = wire.getConnectedElectrodes()
                for electrode in range(numConnected):
                    edges_list.append( ("e" + str(tempList[electrode].getIndex()), w_nodes[-1]))
                
        #print(w_nodes)
        #print(edges_list)
        """
        store numbers of valid wires and edges for generating
        random bipartite graphs later. Note: assumed that graph
        is connected (no isolated nodes), and since only wires that
        connected multiple electrodes were selected, the only way this
        assumption is wrong is if an electrode is not connected to another
        by any wires. High probability that this never happens, so numValidElectrodes
        is just self.num_e
        
        """
        self.numValidWires = len(w_nodes)
        self.numValidEdges = len(edges_list)
        #print("Number of valid wires: ", self.numValidWires)
        #print("Number of valid edges: ", self.numValidEdges)

        # generate bipartite graph
        [graph, _] = generateBipartiteGraph(e_nodes, w_nodes, edges_list)
        return [graph, _]
    
    # generate equivalent bipartite graph
    def generateGraph(self):

        # initialize graph
        g = nx.Graph()

        # add electrode nodes
        for electrode in self.ElectrodeList:
            i = electrode.getIndex()
            x = i % (math.sqrt(self.num_e))
            y = i // (math.sqrt(self.num_e))
            g.add_node(i, pos=(x,y))
            g.add_node(i)
        
        # add wire segments
        for wire in self.Wires:
            # only care if the wire connects 2 or more electrodes
            e_list = wire.getConnectedElectrodes()
            numConnected = len(e_list)
            if numConnected >= 2:
                # outer loop: node 1
                # don't need to worry about last electrode,
                # all possible connections will have been made already
                for e1 in range(numConnected - 1):
                    # inner loop: node 2
                    for e2 in range(e1, numConnected):
                        g.add_edge(e_list[e1].getIndex(), e_list[e2].getIndex())
        
        # draw graph figure
        pos = nx.get_node_attributes(g, 'pos')
        nx.draw(g,pos, with_labels=True)
        #plt.show()
        plt.close()
        
        # return graph object
        return g

        





















        

