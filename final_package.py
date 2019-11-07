from Readfile import *
from Edmonds import *

with open("param") as data:
    line = data.read().split("\n")
    filename = line[3]
    writename = line[6]

nodes = read(filename)
node_name_list = nodes.keys()
(g, root) = create_tree(nodes, node_name_list)
print("finished creating tree.")
result = compute_rdmst(g, "root")

with open(writename, "w+") as write:
    write.write(result.__str__())
