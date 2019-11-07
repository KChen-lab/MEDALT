## Qihan Wang
## qw15
## COMP 182 Homework 3 Problem 3

from collections import *
from copy import *
import pickle

def reverse_digraph_representation(graph):
    """
    Gives a reversed representation of the graph imported

    Arguments:
    graph -- a weighted digraph in standard dictionary representation

    Return:
    A reversed representation of the digraph
    """
    rev_graph = {}
    for node in graph:  # initialize all nodes in the reversed graph
        rev_graph[node] = {}
    for a_node in graph:
        for b_node in graph[a_node]:
            rev_graph[b_node][a_node] = graph[a_node][b_node]
    return rev_graph
# print reverse_digraph_representation(homework3testgraphs.g0)
# print "expected: {0: {}, 1: {0: 2, 4: 2}, 2: {0: 2, 1: 2}, 3: {0: 2, 2: 2}, 4: {2: 2, 3: 2}, 5: {1: 2, 3: 2}}"
# print reverse_digraph_representation(homework3testgraphs.g1)
# print "expected: {0: {}, 1: {0: 20, 4: 4}, 2: {0: 4, 1: 2}, 3: {0: 20, 2: 8}, 4: {2: 20, 3: 4}, 5: {1: 16, 3: 8}}"
# print reverse_digraph_representation(homework3testgraphs.g2)
# print "expected: {0: {}, 1: {0: 5, 2: 2}, 2: {0: 4, 1: 2}}"


def modify_edge_weights(rgraph, root):
    """
    Modifies the edge weights of graph according to Lemma 2

    Arguments:
    rgraph -- a weighted digraph in reversed dictionary representation
    root -- a node in graph.

    Modifies:
    the weight of graph

    Returns:
    No return value
    """
    for node in rgraph:
        if not node == root:
            min = float('inf')
            for in_nbr in rgraph[node]:  # get the minimum edge weight
                if rgraph[node][in_nbr] < min:
                    min = rgraph[node][in_nbr]
            for in_nbr in rgraph[node]:  # modify edge weight
                rgraph[node][in_nbr] -= min
# g = homework3testgraphs.g0
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# print g
# print "expected: {0: {}, 1: {0: 0, 4: 0}, 2: {0: 0, 1: 0}, 3: {0: 0, 2: 0}, 4: {2: 0, 3: 0}, 5: {1: 0, 3: 0}}"
# g = homework3testgraphs.g1
# g = reverse_digraph_representation(g)
# print g
# print "expected: {0: {}, 1: {0: 20, 4: 4}, 2: {0: 4, 1: 2}, 3: {0: 20, 2: 8}, 4: {2: 20, 3: 4}, 5: {1: 16, 3: 8}}"
# g = homework3testgraphs.g2
# g = reverse_digraph_representation(g)
# print g
# print "expected: {0: {}, 1: {0: 5, 2: 2}, 2: {0: 4, 1: 2}}"


def compute_rdst_candidate(rgraph, root):
    """
    Computes an RDST candidate based on Lemma 1 of the weighted digraph.

    Arguments:
    rgraph -- a weighted digraph in reversed dictionary representation
    root -- a node in graph.

    Returns:
    The RDST candidate as a weighted digraph in reverse representation
    """
    candidate = {}
    for node in rgraph:  # initialize all nodes in the reversed graph
        candidate[node] = {}
    for node in rgraph:
        if not node == root:
            min = float('inf')
            for in_nbr in rgraph[node]:  # get the minimum edge weight
                if rgraph[node][in_nbr] < min:
                    min = rgraph[node][in_nbr]
            for in_nbr in rgraph[node]:  # put edge in candidate graph
                if candidate[node] == {}:  # if the node didn't already have an edge
                    if rgraph[node][in_nbr] == min:
                        candidate[node][in_nbr] = min
    return candidate
# g = homework3testgraphs.g0
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# print compute_rdst_candidate(g, 0)
# print "expected: {0: {}, 1: {0: 0}, 2: {0: 0}, 3: {0: 0}, 4: {2: 0}, 5: {1: 0}}"
# g = homework3testgraphs.g1
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# print compute_rdst_candidate(g, 0)
# print "expected: {0: {}, 1: {4: 0}, 2: {1: 0}, 3: {2: 0}, 4: {3: 0}, 5: {3: 0}}"
# g = homework3testgraphs.g2
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# print compute_rdst_candidate(g, 0)
# print "expected: {0: {}, 1: {2: 0}, 2: {1: 0}}"


def compute_cycle(rdst_candidate):
    """
    Computes a cycle in a weighted digraph

    Arguments:
    rdst_candidate -- a RUST candidate as a weighted digraph in reverse representation

    Return:
    A tuple that contains the nodes on the cycle organized according to the cycle
    False if there are no loops in the graph
    """
    node_unvisited = []
    for node in rdst_candidate:  # mark all nodes unvisited
        node_unvisited.append(node)
    while not node_unvisited == []:
        start_node = node_unvisited.pop()  # start DFS on a node
        stack = []  # stack for DFS to run
        trail = []  # keep track of the trail
        stack.append(start_node)
        while not len(stack) == 0:
            node = stack.pop(-1)
            for nbr in rdst_candidate[node]:
                if nbr in trail:
                    return tuple(trail[trail.index(nbr):])  # give the circle
                else:
                    stack.append(nbr)
                    trail.append(nbr)
                    if nbr in node_unvisited:
                        node_unvisited.remove(nbr)
    return False
# g = homework3testgraphs.g0
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# print compute_cycle(rdst)
# print "expected: False"
# g = homework3testgraphs.g1
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# print compute_cycle(rdst)
# print "expected: (3, 2, 1, 4)"
# g = homework3testgraphs.g2
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# print compute_cycle(rdst)
# print "expected: (1, 2)"


def contract_cycle(graph, cycle):
    """
    Contracts the circle found into a node and forms a new graph

    Arguments:
    graph --  a weighted digraph in standard dictionary representation
    cycle -- a tuple that contains the nodes on the cycle organized according to the cycle

    Return:
    contracted_graph -- a contracted digraph with the circle replaced by a node in standard representation
    cstar -- the number of the new node added to replace the cycle
    """
    cstar = max(graph.keys()) + "1"
    contracted_graph = {}
    contracted_graph[cstar] = {}
    for node in graph:  # initialize all nodes in the new graph
        if not node in cycle:
            contracted_graph[node] = {}
    for node in graph:  # add the edges into the new graph
        for nbr in graph[node]:
            if node in cycle:
                if nbr in cycle:
                    pass
                else:  # update the minimum edge weight
                    if nbr in contracted_graph[cstar]:
                        contracted_graph[cstar][nbr] = min(contracted_graph[cstar][nbr], graph[node][nbr])
                    else:
                        contracted_graph[cstar][nbr] = graph[node][nbr]
            else:
                if nbr in cycle:
                    if cstar in contracted_graph[node]:
                        contracted_graph[node][cstar] = min(contracted_graph[node][cstar], graph[node][nbr])
                    else:
                        contracted_graph[node][cstar] = graph[node][nbr]
                else:
                    contracted_graph[node][nbr] = graph[node][nbr]
    return contracted_graph, cstar
# g = homework3testgraphs.g1
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# cycle = compute_cycle(rdst)
# g = reverse_digraph_representation(g)
# print contract_cycle(g, cycle)
# print "expected: ({0: {6: 2}, 5: {}, 6: {5: 0}}, 6)"
# g = homework3testgraphs.g2
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# cycle = compute_cycle(rdst)
# g = reverse_digraph_representation(g)
# print contract_cycle(g, cycle)
# print "expected: ({0: {3: 2}, 3: {}}, 3)"
# g = homework3testgraphs.g3
# g = reverse_digraph_representation(g)
# modify_edge_weights(g, 0)
# rdst = compute_rdst_candidate(g, 0)
# cycle = compute_cycle(rdst)
# g = reverse_digraph_representation(g)
# print contract_cycle(g, cycle)
# print "expected: ({1: {2: 1.1, 4: 0.0, 6: 1.0}, 2: {1: 1.1, 4: 7.9, 6: 1.0}, 4: {1: 8.1, 2: 16.1, 6: 16.0}, 6: {1: 0.0, 2: 0.0, 4: 6.9}}, 6)"


def expand_graph(graph, rdst_candidate, cycle, cstar):
    """
    Expands the contracted digraph to the original form

    Arguments:
    graph -- a weighted digraph in standard representation whose cycle was contracted
    rdst_candidate --  the RDST candidate as a weighted digraph in standard representation
    cycle -- a tuple that contains the nodes on the cycle organized according to the cycle
    cstar -- the number of the new node added to replace the cycle

    Return:
    A weighted digraph in standard representation that results from expanding the cycle in rdst_candidate
    """
    restored_graph = {}
    for node in graph:  # initialize an empty graph
        restored_graph[node] = {}
    for node in rdst_candidate:  # restore edges in rdst_candidate
        for nbr in rdst_candidate[node]:
            if node == cstar:  # find the original node that has an edge to nbr
                min = float('inf')
                for orig in cycle:
                    if nbr in graph[orig]:
                        if graph[orig][nbr] < min:
                            min = graph[orig][nbr]
                            point = orig
                restored_graph[point][nbr] = min
            else:
                if nbr == cstar:  # find the node that origin has an edge to
                    min = float('inf')
                    for orig_nbr in graph[node]:
                        if orig_nbr in cycle:
                            if graph[node][orig_nbr] < min:
                                min = graph[node][orig_nbr]
                                start_pt = orig_nbr  # keep track of the starting point
                    restored_graph[node][start_pt] = min
                else:  # add all other nodes
                    restored_graph[node][nbr] = graph[node][nbr]
    for index in range(len(cycle) - 1):  # restore nodes in cycle
        restored_graph[cycle[index + 1]][cycle[index]] = graph[cycle[index + 1]][cycle[index]]
    restored_graph[cycle[0]][cycle[-1]] = graph[cycle[0]][cycle[-1]]  # add last edge in the cycle
    index = cycle.index(start_pt)
    if index == len(cycle) - 1:
        restored_graph[cycle[0]].pop(cycle[index])
    else:
        restored_graph[cycle[index + 1]].pop(cycle[index])
    return restored_graph
# (graph, rdst_candidate, cycle, cstar) = \
#     ({1: {8: 0.0, 2: 1.1, 3: 1.0, 5: 1.1, 14: 0.0}, 2: {8: 8.0, 1: 2.1, 3: 1.0, 5: 1.0, 14: 7.9},
#       3: {8: 7.0, 1: 1.0, 2: 0.0, 5: 0.0, 14: 6.9}, 5: {8: 7.0, 1: 1.1, 2: 0.0, 3: 0.0, 14: 6.9},
#       8: {1: 6.1, 2: 13.1, 3: 13.1, 5: 13.1, 14: 6.0}, 14: {8: 9.0, 1: 9.1, 2: 16.1, 3: 16.0, 5: 16.0}},
#      {8: {}, 1: {8: 0.0, 14: 0.0, 15: 0.0}, 2: {}, 14: {}, 15: {2: 0.0}}, (3, 5), 15)
# print expand_graph(graph, rdst_candidate, cycle, cstar)
# print "expected: {1: {8: 0.0, 3: 1.0, 14: 0.0}, 2: {}, 3: {2: 0.0, 5: 0.0}, 5: {}, 8: {}, 14: {}}"
# (graph, rdst_candidate, cycle, cstar) = \
#     ({1: {2: 1.1, 3: 1.0, 5: 1.1, 8: 0.0, 9: 4.9, 13: 3.0},
#       2: {1: 2.1, 3: 1.0, 5: 1.0, 8: 8.0, 9: 13.000000000000002, 13: 10.9},
#       3: {1: 1.0, 2: 0.0, 5: 0.0, 8: 7.0, 9: 12.000000000000002, 13: 9.9},
#       5: {1: 1.1, 2: 0.0, 3: 0.0, 8: 7.0, 9: 12.000000000000002, 13: 9.9},
#       8: {1: 6.1, 2: 13.1, 3: 13.1, 5: 13.1, 9: 11.000000000000002, 13: 9.0},
#       9: {1: 11.1, 2: 18.1, 3: 18.1, 5: 18.1, 8: 11.000000000000002, 13: 0.0},
#       13: {1: 9.1, 2: 16.1, 3: 16.0, 5: 16.0, 8: 9.0, 9: 0.0}},
#      {1: {8: 0.0, 3: 1.0, 14: 0.0}, 2: {}, 3: {2: 0.0, 5: 0.0}, 5: {}, 8: {}, 14: {}}, (9, 13), 14)
# print expand_graph(graph, rdst_candidate, cycle, cstar)
# print "expected: {1: {8: 0.0, 3: 1.0, 13: 3.0}, 2: {}, 3: {2: 0.0, 5: 0.0}, 5: {}, 8: {}, 9: {}, 13: {9: 0.0}}"
# (graph, rdst_candidate, cycle, cstar) = \
#     ({1: {2: 1.1, 3: 1.0, 4: 4.0, 5: 1.1, 8: 0.0, 9: 4.9, 12: 5.1},
#       2: {1: 2.1, 3: 1.0, 4: 11.9, 5: 1.0, 8: 8.0, 9: 13.000000000000002, 12: 13.0},
#       3: {1: 1.0, 2: 0.0, 4: 10.9, 5: 0.0, 8: 7.0, 9: 12.000000000000002, 12: 12.0},
#       4: {1: 9.1, 2: 16.1, 3: 16.0, 5: 16.0, 8: 9.0, 9: 0.0, 12: 0.0},
#       5: {1: 1.1, 2: 0.0, 3: 0.0, 4: 10.9, 8: 7.0, 9: 12.000000000000002, 12: 12.0},
#       8: {1: 6.1, 2: 13.1, 3: 13.1, 4: 10.0, 5: 13.1, 9: 11.000000000000002, 12: 11.0},
#       9: {1: 11.1, 2: 18.1, 3: 18.1, 4: 1.0, 5: 18.1, 8: 11.000000000000002, 12: 2.0},
#       12: {1: 10.1, 2: 17.1, 3: 17.0, 4: 0.0, 5: 17.0, 8: 9.9, 9: 0.9000000000000004}},
#      {1: {8: 0.0, 3: 1.0, 13: 3.0}, 2: {}, 3: {2: 0.0, 5: 0.0}, 5: {}, 8: {}, 9: {}, 13: {9: 0.0}}, (4, 12), 13)
# print expand_graph(graph, rdst_candidate, cycle, cstar)
# print "expected: " \
#       "{1: {8: 0.0, 3: 1.0, 4: 4.0}, 2: {}, 3: {2: 0.0, 5: 0.0}, 4: {9: 0.0, 12: 0.0}, 5: {}, 8: {}, 9: {}, 12: {}}"


def bfs(graph, startnode):
    """
        Perform a breadth-first search on digraph graph starting at node startnode.

        Arguments:
        graph -- directed graph
        startnode - node in graph to start the search from

        Returns:
        The distances from startnode to each node
    """
    dist = {}

    # Initialize distances
    for node in graph:
        dist[node] = float('inf')
    dist[startnode] = 0

    # Initialize search queue
    queue = deque([startnode])

    # Loop until all connected nodes have been explored
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def compute_rdmst(graph, root):
    """
        This function checks if:
        (1) root is a node in digraph graph, and
        (2) every node, other than root, is reachable from root
        If both conditions are satisfied, it calls compute_rdmst_helper
        on (graph, root).

        Since compute_rdmst_helper modifies the edge weights as it computes,
        this function reassigns the original weights to the RDMST.

        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node id.

        Returns:
        An RDMST of graph rooted at r and its weight, if one exists;
        otherwise, nothing.
    """

    if root not in graph:
        print "The root node does not exist"
        return

    distances = bfs(graph, root)
    for node in graph:
        if distances[node] == float('inf'):
            print "The root does not reach every other node in the graph"
            return

    rdmst = compute_rdmst_helper(graph, root, 1)

    # reassign the original edge weights to the RDMST and computes the total
    # weight of the RDMST
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = graph[node][nbr]
            rdmst_weight += rdmst[node][nbr]

    return (rdmst,rdmst_weight)


def compute_rdmst_helper(graph, root, index):
    """
        Computes the RDMST of a weighted digraph rooted at node root.
        It is assumed that:
        (1) root is a node in graph, and
        (2) every other node in graph is reachable from root.

        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node in graph.

        Returns:
        An RDMST of graph rooted at root. The weights of the RDMST
        do not have to be the original weights.
        """

    # reverse the representation of graph
    rgraph = reverse_digraph_representation(graph)

    # Step 1 of the algorithm
    modify_edge_weights(rgraph, root)

    # Step 2 of the algorithm
    rdst_candidate = compute_rdst_candidate(rgraph, root)

    # compute a cycle in rdst_candidate
    cycle = compute_cycle(rdst_candidate)

    # Step 3 of the algorithm
    if not cycle:
        return reverse_digraph_representation(rdst_candidate)
    else:
        # Step 4 of the algorithm

        g_copy = deepcopy(rgraph)
        g_copy = reverse_digraph_representation(g_copy)

        # Step 4(a) of the algorithm
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        #cstar = max(contracted_g.keys())

        with open(str(index) + '.txt', 'wb') as cache:
            pickle.dump(rgraph, cache)

        del graph, rgraph, rdst_candidate, g_copy

        print("contracting cycle " + str(index))

        # Step 4(b) of the algorithm
        new_rdst_candidate = compute_rdmst_helper(contracted_g, root, index + 1)

        with open(str(index) + '.txt', 'rb') as cache:
            rgraph = pickle.loads(cache.read())

        print("expanding cycle " + str(index))

        # Step 4(c) of the algorithm
        rdmst = expand_graph(reverse_digraph_representation(rgraph), new_rdst_candidate, cycle, cstar)

        return rdmst
