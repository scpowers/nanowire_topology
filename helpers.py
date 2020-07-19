import networkx as nx
import matplotlib.pyplot as plt

# helpers.py: defines helper functions for use in other classes. 
# mostly comprised of methods for graph analysis

# generate equivalent bipartite graph from geometric model
def generateBipartiteGraph(e_node_set, w_node_set, edge_set):
    
    # initialize graph
    graph = nx.Graph()

    # add nodes to the graph
    graph.add_nodes_from(e_node_set, bipartite=0) # electrode node set
    graph.add_nodes_from(w_node_set, bipartite=1) # wire node set

    # add edges
    graph.add_edges_from(edge_set)

    # create dictionary for bipartite graph layout
    pos = nx.bipartite_layout(graph, nx.bipartite.sets(graph)[0])

    # assign colors to the two node sets
    color_dict = {0:'r', 1:'b'}
    color_list = [color_dict[i[1]] for i in graph.nodes.data('bipartite')]

    # draw graph
    nx.draw(graph, pos, with_labels=True, node_color=color_list)
    #plt.show()
    plt.close()

    # return graph and set of electrode nodes (for analysis)
    return [graph, nx.bipartite.sets(graph)[0]]

# generate random bipartite graph (needed for analysis)
def generateRandomBipartiteGraph(n, m, k):
    # create new NMK graph:
    # n = number of nodes in first bipartite set (numElectrodes)
    # m = number of nodes in second bipartite set (numWires)
    # k = number of edges (numValidEdges)
    graph = nx.bipartite.gnmk_random_graph(n, m, k)

    while (not nx.is_connected(graph)):
        graph = nx.bipartite.gnmk_random_graph(n, m, k)

    # isolate top nodes for pos dictionary definition
    top = [node for node in graph.nodes() if graph.nodes[node]['bipartite']==0]

    # create dictionary for bipartite graph layout
    pos = nx.bipartite_layout(graph, top)

    # assign colors to the two node sets
    color_dict = {0:'r', 1:'b'}
    color_list = [color_dict[i[1]] for i in graph.nodes.data('bipartite')]

    # draw graph
    nx.draw(graph, pos, node_color=color_list)
    #plt.show()
    plt.close()
    
    # return graph and set of equivalent electrode nodes for analysis
    return [graph, top]
    
# analyze bipartite graph
def analyzeBipartiteGraph(graph, e_nodes):

    # average shortest path length
    L_r = nx.average_shortest_path_length(graph)
    # square clustering coefficient of electrodes
    C_r_dict = nx.square_clustering(graph, e_nodes)
    C_r_list = []
    for key in C_r_dict:
        C_r_list.append(C_r_dict[key])
    C_r = sum(C_r_list) / len(C_r_list)

    # return L_r and C_r for later analysis
    return [L_r, C_r]
    
# compute small world coefficient, given geometric model bipartite
# graph and a random bipartite graph with the same number of nodes
# in each set and the same number of edges
def computeSmallWorldCoefficient(L, C, L_r, C_r):
    return ( C / C_r ) / ( L / L_r ) 


# plot small world coefficient versus square root of electrode count
def plotSmallWorldSweep(sigma_list, e_count_list):
    plt.plot(e_count_list, sigma_list, '-o')
    plt.xlabel(r'$\sqrt{n_e}$')
    plt.ylabel(r'$\sigma$')
    plt.xticks([0, 10, 20, 30, 40, 50, 60])
    plt.savefig('sweep.png')
