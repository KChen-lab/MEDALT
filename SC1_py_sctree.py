# required packages
import os,sys
import numpy as np
import pandas as pd
from datetime import datetime as dt_
from SP1_SCT_UTIL import *

def main():
    ### Initialing variables #####################################################################
    ##############################################################################################
    sys.setrecursionlimit(20000)
    usage = "usage: python %prog <-P path> <-I input> <-D datatype>"
    description = """Input integer copy number profile, diploid conversion should be done prior to this run. 
                     Columns correspond to cell IDs and Rows correspond to genomic regions (scDNA) or genes (scRNA)."""
    op = OptionParser(version="%prog 1.0", description=description, usage=usage, add_help_option=False)
    op.add_option("-h", "--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-P", "--Path",dest="Path",type="str",
                  help="Path to the script")
    op.add_option("-I", "--Input",dest="Input",type="str",
                  help="Input file")
    op.add_option("-G", "--Genome",dest="Genome",type="str",
                  help="Genome version hg19 or hg38")
    op.add_option("-O", "--Output",dest="Output",type="str",
                  help="Output path. Default: input name + runtime")
    op.add_option("-D", "--Datatype",dest="Datatype",type="str",
                  help="The type of input data. Either D (DNA-seq) or R (RNA-seq).")
    op.add_option("-W", "--Windows",dest="Windows",type="str",
                  help="the number of genes to merge when the input data type is R. Default 30.")
    op.add_option("-R", "--Permutation",dest="Permutation",type="str",
                  help="""Whether reconstructed permuted tree (T) or not (F). 
                          If not, permuted copy number profile will be used to perform LSA. 
                          Default value is F due to time cost.""")

    (options,args) = op.parse_args()
    # check input parameters. Package path, input file, data type and genome version are required.
    if not options.Path or not options.Input or not options.Datatype or not options.Genome:
        op.print_help() ; sys.exit(1)
    
    # get the input parameters
    PCKAGE_PATH = options.Path
    IN_CNV_PATH = options.Input
    IN_CNV_FILE = IN_CNV_PATH.split("/")[-1][:-4]
    NUCLEC_ACID = options.Datatype
    REF__GENOME = options.Genome
    GENE_BIN_SZ = options.Windows if options.Windows else "30"
    rt = dt_.now() ; rt = f"{rt.year}_{rt.month:0>2}{rt.day:0>2}_{rt.hour:0>2}{rt.minute:0>2}"
    OUTPUT_PATH = options.Output if options.Output else f"{IN_CNV_PATH[:-4]}_{rt}"
    if not options.Permutation: permutation = "F"
    else: permutation = options.Permutation

    # get the input parameters
    PCKAGE_PATH = getPath(PCKAGE_PATH).replace("//", "/")
    IN_CNV_PATH = getPath(IN_CNV_PATH).replace("//", "/")
    OUTPUT_PATH = getPath(OUTPUT_PATH).replace("//", "/")
    GENPOS_PATH = f"{PCKAGE_PATH}/genomes/gencode_v{REF__GENOME[-2:]}_gene_pos.txt"
    BANDBD_PATH = f"{PCKAGE_PATH}/genomes/{REF__GENOME}.band.bed"

    # supplementary files paths
    DE_DUP_PATH = f"1_{IN_CNV_FILE}_dedup.csv"
    DUPREF_PATH = f"1_{IN_CNV_FILE}_dup_ref.csv"
    SEGCNV_PATH = f"2_{IN_CNV_FILE}_bin_{GENE_BIN_SZ}.csv"
    SCTREE_PATH = f"3_CNV.tree.txt"
    os.system("mkdir -p " + OUTPUT_PATH)

    
    ### Data reforming ###########################################################################
    ### change the input copy number profile to the matrix format used to infer Tree           ###
    ### if the data type = R, estimate integer copy number by averaging adjacent genes.        ###
    ### Default number of genes is 30.                                                         ###
    ##############################################################################################
    # import cell CNV data
    print("\n#####################################################")
    print("### now running SC1_py_sctree.py                  ###")
    print("#####################################################\n")
    os.chdir(OUTPUT_PATH)
    print(f"All intermediate files and output files will be stored in {OUTPUT_PATH}")
    print("reading data")
    df_ori = pd.read_csv(IN_CNV_PATH, sep="\t", index_col=0)
    df_ori.columns = df_ori.columns.str.replace("[\ \-\.]", "_", regex=True)

    # deduplication
    cel_dup = df_ori.T.duplicated(keep="first")
    cel_uni = cel_dup[~cel_dup].index.tolist()  # unique cell list
    cel_dup = cel_dup[ cel_dup].index.tolist()  # duplicated cell list
    # find duplication relationship
    print("running deduplication")
    dup_relationship = {"par_cell":[], "dup_cell":[]}
    for dup_cell in cel_dup:
        for uni_cell in cel_uni:
            if df_ori[uni_cell].equals(df_ori[dup_cell]):  # find cells of origin
                dup_relationship["par_cell"].append(uni_cell)
                dup_relationship["dup_cell"].append(dup_cell)
                break
    cell_dup_ref = pd.DataFrame(dup_relationship)
    if cell_dup_ref.shape[0]>0:
        # save duplication information
        df_ded = df_ori.loc[:,cel_uni]
        #df_ded = df_ded.iloc[:,np.random.choice(range(df_ded.shape[1]), size=min(df_ded.shape[1], 100), replace=False)]
        print(f"""{df_ded.shape[1]}/{df_ori.shape[1]} cells remained after deduplication.
                        duplication relationship is stored in {DUPREF_PATH} and 
                        dedup-data saved in {DE_DUP_PATH}, which will be used in the following analysis""")
        df_ded.to_csv(DE_DUP_PATH, sep="\t")
        cell_dup_ref.to_csv(DUPREF_PATH, sep="\t")
    else:
        print("no duplicating cells found, good!")
        DE_DUP_PATH = IN_CNV_PATH
        df_ded = df_ori
    
    print(f"converting CN profile to segmental CN level:\n")
    command = f"Rscript {PCKAGE_PATH}/SC2_RR_dataTransfer.R {OUTPUT_PATH} {DE_DUP_PATH} {NUCLEC_ACID} {SEGCNV_PATH}"                                # input DNA
    if NUCLEC_ACID == "R": command = f"{command} {GENPOS_PATH} {GENE_BIN_SZ}"   # input RNA
    os.system(command) ; print(command)

    ### MEDLAT tree inference ####################################################################
    ### reformed segmental data will be used to infer MEDLAT single cell tree                  ###
    ### where a diploid cell will be used as the root (imputed if not existing)                ###
    ##############################################################################################
    nodes, root = read_CNV(SEGCNV_PATH)
    node_list  = nodes.keys()

    #calculation of MED distance
    print("\n#####################################################")
    print("### going back to SC1_py_sctree.py                ###")
    print("#####################################################\n")
    print("initializing tree")
    tree_dict = create_tree(nodes, node_list, root, df_cor=None, len_threshold=30)  
    # set df_cor to None and leave proximity to True if no spatial coordinate information is provided
    # this will automatically calculate pairwise MED instead of only connecting cells within close proximity
    
    print("computing rdmst")
    tree = compute_rdmst(tree_dict, root)[0]
    with open(SCTREE_PATH,'w') as write:
        write.write("\t".join(["stt", "end", "len"])+"\n") # header line
        for in_node in tree.keys():
            for ot_node in tree[in_node].keys():
                write.write("\t".join([in_node, ot_node, str(tree[in_node][ot_node])])+"\n")
    print(f"MEDALT inferrence finish, weighted tree saved to:{SCTREE_PATH}")

    #Permutation process for lineage speciation analysis (LSA)
    ### LSA test #################################################################################
    ### lineage speciation analysis can be done with/without permutation.                      ###
    ### with permutation (T), both copy number profiles and single cells labels are shuffled   ###
    ### to repeat above inference process and be used to estimate LSA result significance.     ###
    ### otherwise (permu==F), only shuffle single cell labels to estimate significance.        ###
    ##############################################################################################
    if permutation == "T":
        PERMUT_PATH = OUTPUT_PATH + "/permutation"
        print("Reconstructing tree based on permutation data.")
        print("This will take a long time! Please have some coffee.")

        #permute copy number profile
        step_3_cmd = f"Rscript {PCKAGE_PATH}/SC3_RR_Permutation.R {SCTREE_PATH} "
        step_3_cmd = step_3_cmd + f"{IN_CNV_PATH} {NUCLEC_ACID} {OUTPUT_PATH}/permutation"
        if datatype == "D":
            os.system(step_3_cmd)
        elif datatype == "R":
            os.system(step_3_cmd + " "+ GENE_BIN_SZ + " " + REF__GENOME)

        #Infer permutation tree
        for j in range(1,101):
            permutefile=permutationPath+"/permute."+str(j)+".CNV.txt"
            (nodes,root) = read(permutefile)
            node_name_list = nodes.keys()
            g = create_tree(nodes, node_name_list,root)
            result = compute_rdmst(g, root)
            permuteTree=permutefile+".celltree.txt"
            write=open(permuteTree,'w')
            tree=result[0]
            out1="from"+"\t"+"to"+"\t"+"dist"
            write.write(out1)
            for ele in tree.keys():
                out=ele
                for value in tree[ele].keys():
                    out1=out+"\t"+value+"\t"+str(tree[ele][value])
                    print >> write,out1
            write.close()
        print("Pemutation tree finish.")
    elif permutation == "F": PERMUT_PATH = ""

    #Identifying CNAs associated with cellular lineage expansion.
    print("Performing LSA.")
    cmd = f"Rscript {PCKAGE_PATH}/SC4_RR_LSA_tree.R {PCKAGE_PATH} {IN_CNV_PATH} {OUTPUT_PATH} {SCTREE_PATH} {SEGCNV_PATH} {NUCLEC_ACID} {REF__GENOME} {PERMUT_PATH}"
    print(cmd)
    os.system(cmd)
    print("All is done!")


if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
