##############################################################################################################################
# Returns a dictionary mapping node names to list of list of integers representing list of copy number list
##############################################################################################################################
def read_CNV(in_seg_path):
    df = pd.read_csv(in_seg_path, sep="\t")  # read data
    
    chr_scan = [f"chr{i}" for i in range(1,25)] + ["chrX", "chrY"]  # def candidate chrs
    chr_exis = df.columns.str.replace("_.*$", "", regex=True)       # find existing chrs
    chr_locs={}
    for chr_i in chr_scan:
        chr_segs = sum(chr_exis==chr_i)         # chr seg number
        if chr_segs>0: chr_locs[chr_i]=chr_segs # attach it to the dict
    cell_lst = list(df.index)                   # cell ids
    node_dic = {i:[]    for i in cell_lst}      # to store node info
    cnvs_dic = {i:set() for i in cell_lst}      # to keep track of CNs
    for chr_i in chr_locs:
        seg_i = chr_exis==chr_i
        for cel_j in cell_lst: 
            node_dic[cel_j] = node_dic[cel_j] +[list(df.loc[cel_j, seg_i])]# concat chr bin CNs and append
            cnvs_dic[cel_j] = cnvs_dic[cel_j] |  set(df.loc[cel_j, seg_i]) # keep track of CNV in each cell
    
    root = False
    for cel_j in cell_lst:
        if cnvs_dic[cel_j]=={2}: root = cel_j
    if not root:
        print(f"No diploid found, inputating a root cell.")
        root='root'
        node_dic[root] = [[2]*chr_locs[i] for i in chr_locs.keys()]
    return node_dic, root

##############################################################################################################################

##############################################################################################################################
def getPath(path):
    abs_path = os.getcwd().split("/")
    for i,ele in enumerate(path.split("/")):
        if ele == ".":
            pass
        elif ele == "..":
            abs_path = abs_path[:-1]
        elif (ele == "") & (i == 0):
            abs_path = path.split("/")
            break
        else: 
            abs_path = abs_path + [ele]
    return "/".join(abs_path)


##############################################################################################################################

##############################################################################################################################
#class Tree:
#    def __init__(self, name, ot_edge, in_edge):
#        self.name = name  # name of node
#        self.ot_edge = ot_edge  # dictionary of node name - distance pair
#        self.in_edge = in_edge  # dictionary of node name - distance pair
#
#    def get_out_degree(self):
#        return len(self.ot_edge)
#
#    def get_in_degree(self):
#        return len(self.in_edge)
#    
###############################################################################################################################

##############################################################################################################################
def create_tree(nodes, node_list, root, proximity=True, len_threshold=30, df_cor=None):
    print("{:4d} cells to run.".format((len(node_list))), end="")
    blk = 5 ; sep = max(5,len(node_list)//blk)
    tp0 = dt_.now() ; tpt = dt_.now() # for timing
    tree_node_dict = {}
    for idx, A_node in enumerate(node_list, 1):
        pct = min(1, idx/len(node_list))  ; ext = 1/pct-1                   # used to report runtime
        nwl = "\n" if (idx-1)%sep==0 else "\r"                              # used to report runtime
        dtt = dt_.now()-tp0 ; rmt = str(dtt*ext)[2:-5] ; dtt=str(dtt)[2:-5] # used to report runtime
        print(f"""{nwl} iter {idx:4d}, {pct*100:5.1f}% , {
                ""}elapse {dtt}, expected in {rmt}""", end="")              # used to report runtime
            
        out_edge = {} ; candidates = list(set(node_list)^set([A_node]))
        for B_node in candidates:
            if df_cor is not None:
                physical_dist = ((df_cor.loc[A_node, "coor_x"] - df_cor.loc[B_node, "coor_x"])**2+(df_cor.loc[A_node, "coor_y"] - df_cor.loc[B_node, "coor_y"])**2)
                proximity = True if physical_dist<len_threshold**2 else False
            if proximity:
                if B_node in tree_node_dict.keys():
                    if A_node in tree_node_dict[B_node].keys():
                        out_edge[B_node] = tree_node_dict[B_node][A_node]
                else:   out_edge[B_node] = dist(nodes[A_node], nodes[B_node])
                #out_edge[B_node] = dist(nodes[A_node], nodes[B_node])
            else: 
                pass # out_edge[B_node] = 9999
        tree_node_dict[A_node] = out_edge
    print("\ntotal tree initiation time: {}".format(dt_.now()-tp0))
    # final sanity check of the tree
    if root not in tree_node_dict: print(f"the root node: {root}, not found in the graph")
    distances =  graph_rank_dist(tree_node_dict, root)
    for node in tree_node_dict:
        if distances[node] == float('inf'): print(f"node {node} not connected to the root, please check it.")
    # no problem!
    return tree_node_dict

##############################################################################################################################

##############################################################################################################################
def matrixbuilder(node):
    matrix = []
    for node1 in node:
        temp = []
        for node2 in node:
            temp.append(dist(node[node1], node[node2]))
        matrix.append(temp)
    return node.keys(), matrix



##############################################################################################################################

##############################################################################################################################
def dist(node1, node2, d=0):
    for i in range(0, len(node1)):
        d = d + distcalc(node1[i], node2[i])
    return d



##############################################################################################################################

##############################################################################################################################
def distcalc(node1, node2):
    assert len(node1) == len(node2)
    d = 0
    diff = copy.deepcopy(node1) ##########################
    for i in range(0, len(node2)):
        diff[i] = diff[i] - node2[i]
    while diff:
        if   diff[0] == 0: diff.pop(0)
        elif diff[0] >  0:
            for i in range(0, len(diff)):
                if diff[i] > 0: diff[i] = diff[i] - 1
                else: break
            d = d+1
        elif diff[0] <  0:
            for i in range(0, len(diff)):
                if diff[i] < 0: diff[i] = diff[i] + 1
                else: break
            d = d+1
    return d



##############################################################################################################################

##############################################################################################################################
def graph_rank_dist(g, startnode):
    dist = {}

    for node in g:
        dist[node] = float('inf')
    dist[startnode] = 0

    queue = deque([startnode])

    while queue:
        node = queue.popleft()
        for nbr in g[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist



##############################################################################################################################

##############################################################################################################################
# inner psutil function
def pc_mem():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return round(mem_info.rss/1024/1024/1024, 2)

##############################################################################################################################

##############################################################################################################################
def compute_rdmst(g, root):
    st_mem = pc_mem() ; t0 = dt_.now() ; ts = dt_.now()
    print("Inferring single cell evolution tree:", end="")
    rdmst = rdmst_recursor(g, root, 0, t0, ts, st_mem)
    print(f"\rtotal time elapse {dt_.now()-t0}{' '*20}")
    
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = g[node][nbr]
            rdmst_weight += rdmst[node][nbr]
    return (rdmst, rdmst_weight)


##############################################################################################################################

##############################################################################################################################
def rdmst_recursor(g__inputed, root_node, recurr_idx, t0, ts, st_mem):
    recurr_idx = recurr_idx + 1         # record how many recursions have taken place
    if recurr_idx%5==0: 
        nwl = "\n" if (recurr_idx-5)%100==0 else "\r"
        print(f"""{nwl}{recurr_idx:4d} loops shrunk, avg_span: {str(dt_.now()-ts)[3:-5]
              } elapse: {str(dt_.now()-t0)[:-5]}, mem: {pc_mem()-st_mem:3.2f} GB.""", end="")
        ts = dt_.now()
    g_reversed = reverse__graph(g__inputed)            # reverse graph direction
    g_rev_cont = contract_graph(g_reversed, root_node) # contract each node by min_length
    g_rdst_min = get_rdst_graph(g_rev_cont, root_node) # get the 0_length nodes 
    loop_in_gr = find_graphloop(g_rdst_min)            # find loops in the graph
    
    if not loop_in_gr: 
        print(f"\nall {recurr_idx:4d} loops contracted, now recovering tree:")
        return reverse__graph(g_rdst_min)
    else:   
        g_contract = reverse__graph(g_rev_cont)
        g___pruned, loop_nodes = prune____graph(g_contract, loop_in_gr)
        del g_reversed, g_rev_cont, g__inputed, g_rdst_min
        g_new_rdst = rdmst_recursor(g___pruned, root_node, recurr_idx, t0, ts, st_mem)  # recursion at here
        g_expanded = expand___graph(g_contract, g_new_rdst, loop_in_gr, loop_nodes)
        print(f"\r{recurr_idx:4d} layer remaining.", end="")
        return g_expanded


##############################################################################################################################

##############################################################################################################################
def reverse__graph(g):
    rev_graph = {}
    for node in g:
        rev_graph[node] = {}
    for a_node in g:
        for b_node in g[a_node]:
            rev_graph[b_node][a_node] = g[a_node][b_node]
    return rev_graph


##############################################################################################################################

##############################################################################################################################
def contract_graph(in_graph, root):
    rg = in_graph.copy()
    for node in rg:
        if not node == root:
            minimum = min(rg[node].values())
            for in_nbr in rg[node]:
                rg[node][in_nbr] -= minimum
    return rg


##############################################################################################################################

##############################################################################################################################
def get_rdst_graph(rg, root):
    candidate = {}
    for node in rg:
        candidate[node] = {}
    for node in rg:
        if not node == root:
            minimum = min(rg[node].values())
            for in_nbr in rg[node]:
                if candidate[node] == {}:  
                    if rg[node][in_nbr] == minimum:
                        candidate[node][in_nbr] = minimum
    return candidate


##############################################################################################################################

##############################################################################################################################
def find_graphloop(rdst_candidate):
    node_unvisited = []
    for node in rdst_candidate:  
        node_unvisited.append(node)
    while node_unvisited != []:
        start_node = node_unvisited.pop()  
        stack = []  
        trail = []  
        stack.append(start_node)
        while len(stack) != 0:
            node = stack.pop(-1)
            for nbr in rdst_candidate[node]:
                if nbr in trail:
                    return tuple(trail[trail.index(nbr):]) 
                else:
                    stack.append(nbr)
                    trail.append(nbr)
                    if nbr in node_unvisited:
                        node_unvisited.remove(nbr)
    return False


##############################################################################################################################

##############################################################################################################################
def prune____graph(g, cycle):
    loop_nodes = max(g.keys()) + "1"
    contracted_graph = {}
    contracted_graph[loop_nodes] = {}
    for node in g:
        if not node in cycle:
            contracted_graph[node] = {}
    for node in g:
        for nbr in g[node]:
            if node in cycle:
                if nbr in cycle: pass   ############## node in,  nbr in  cycle ###########################################################
                else:    ############################# node in,  nbr out cycle ###########################################################
                    if nbr in contracted_graph[loop_nodes]:
                        contracted_graph[loop_nodes][nbr] = min(contracted_graph[loop_nodes][nbr], g[node][nbr])
                    else:
                        contracted_graph[loop_nodes][nbr] = g[node][nbr]
            else:
                if nbr in cycle: ##################### node out, nbr in  cycle ###########################################################
                    if loop_nodes in contracted_graph[node]:
                        contracted_graph[node][loop_nodes] = min(contracted_graph[node][loop_nodes], g[node][nbr])
                    else:
                        contracted_graph[node][loop_nodes] = g[node][nbr]
                else: ################################ node out, nbr out cycle ###########################################################
                    contracted_graph[node][nbr] = g[node][nbr]
    return contracted_graph, loop_nodes


##############################################################################################################################

##############################################################################################################################
def expand___graph(g_contract, g_new_rdst, loop_node, loop_repr):
    restored_graph = {}
    for node in g_contract:
        restored_graph[node] = {}
    for node in g_new_rdst:  
        for nbr in g_new_rdst[node]:
            if node == loop_repr:  
                minimum = float('inf')
                for orig in loop_node:
                    if nbr in g_contract[orig]:
                        if g_contract[orig][nbr] < minimum:
                            minimum = g_contract[orig][nbr]
                            point = orig
                restored_graph[point][nbr] = minimum
            else:
                if nbr == loop_repr:  
                    minimum = float('inf')
                    for orig_nbr in g_contract[node]:
                        if orig_nbr in loop_node:
                            if g_contract[node][orig_nbr] < minimum:
                                minimum = g_contract[node][orig_nbr]
                                start_pt = orig_nbr  
                    restored_graph[node][start_pt] = minimum
                else: 
                    restored_graph[node][nbr] = g_contract[node][nbr]
    for index in range(len(loop_node) - 1): 
        restored_graph[loop_node[index + 1]][loop_node[index]] = g_contract[loop_node[index + 1]][loop_node[index]]
    restored_graph[loop_node[0]][loop_node[-1]] = g_contract[loop_node[0]][loop_node[-1]]
    index = loop_node.index(start_pt)
    if index == len(loop_node) - 1:
        restored_graph[loop_node[0]].pop(loop_node[index])
    else:
        restored_graph[loop_node[index + 1]].pop(loop_node[index])
    return restored_graph


##############################################################################################################################

##############################################################################################################################

from datetime import datetime as dt_
from optparse import OptionParser
from collections import *
import psutil

import copy
import pickle
import os,sys
import subprocess
import numpy as np
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
sys.setrecursionlimit(20000)