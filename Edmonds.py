from Readfile import *
from ComputeDistance import *
from copy import *
from mdmst import *

class Tree:
    def __init__(self, name, out_edge, in_edge):
        self.name = name  # name of node
        self.out_edge = out_edge  # dictionary of node name - distance pair
        self.in_edge = in_edge  # dictionary of node name - distance pair

    def get_out_degree(self):
        return len(self.out_edge)

    def get_in_degree(self):
        return len(self.in_edge)


def create_tree(nodes, node_name_list):
    num_calc = 0
    tree_node_dict = {}
    for node in node_name_list:
        temp_out_edge = {}
        for other_node in node_name_list:
            if not node == other_node:
                temp_out_edge[other_node] = dist(nodes[node], nodes[other_node])
                num_calc += 1
                if (num_calc % 10000 == 0):
                    print("did " + str(num_calc) + " computation.")
        if len(temp_out_edge) < 1000:
            tree_node_dict[node] = temp_out_edge
        else:
            tree_node_dict[node] = sorted(temp_out_edge)[:100]
    return tree_node_dict
