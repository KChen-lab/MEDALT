from optparse import OptionParser 
from Readfile import *
from Edmonds import *


def main():
    usage = "usage: python %prog <-I input> <-O output>"
    description = "Input integer copy number profile. Columns correspond to chromosomal position. Rows correspond to cell or clonal ID"
    op = OptionParser(version="%prog 1.0",description=description,usage=usage,add_help_option=False)
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-I","--Input",dest="Input",type="str",
                  help="Input file")
    op.add_option("-O","--Output",dest="Output",type="str",
                  help="Output file")
    (options,args) = op.parse_args()
    filename=options.Input
    writename=options.Output
    (nodes,root) = read(filename)
    node_name_list = nodes.keys()
    g = create_tree(nodes, node_name_list,root)
    result = compute_rdmst(g, root)
    write=open(writename,'w')
    tree=result[0]
    for ele in tree.keys():
        out=ele
        for value in tree[ele].keys():
            out1=out+"\t"+value+"\t"+str(tree[ele][value])
            print >> write,out1
    write.close()


if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
