from Readfile import *
from ComputeDistance import *
from Kruskal import *
from ete3 import PhyloTree


def main():
    node = read()
    keys, matrix = matrixbuilder(node)
    g = Graph(len(keys))
    for i in range(0, len(keys)):
        for j in range(0, len(keys)):
            g.addEdge(i, j, matrix[i][j])
    g.KruskalMST()
    print keys


def matrixtograph():
    return 0

main()