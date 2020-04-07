Description
===========
This package is to infer Minimal Event Distance Aneuploidy Lineage Tree (MEDALT) using integer copy number profile from single cell sequencing technology. It will infer:
* An rooted directed minimal spanning tree (RDMST) to represent aneuploidy evolution of tumor cells
* The focal and broad copy number alterations associated with lineage expansion


System requirements and dependency
==================================
This package runs on Python 2.7.

It also requires R/3.5
to run and has dependency on the R packages:

	igraph, HelloRanges and DescTools.



Installation
============
Please download and copy the distribution to your specific location. If you are cloning from github, ensure that you have git-lfs installed.

For example, if the downloaded distribuition is MEDALT.tar.gz.
	Type 'tar zxvf MEDALT.tar.gz'

Then, run scTree.py in the resulting folder.

Usage
=====
```
Options:
  --version             show program's version number and exit
  -h, --help            Show this help message and exit.
  -P PATH, --PATH=PATH
                        the path of MEDALT package
  -I INPUT, --Input=INPUT
                        the input file is integer copy number profile estimated from scDNA-seq or scRNA-seq
  -O OUTPATH, --outpath=OUTPATH
                        the output path.
  -D DATATYPE, --DATATYPE=DATATYPE     the input file type either D (scDNA-seq) or R (scRNA-seq)
  -W WINDOWS, --WINDOWS=WINDOWS
                        The size of smoothing windows if your inputfile is from scRNA-seq.
                        The value is the number of genes which will be merge. Default value is 30.
	-W WINDOWS, --WINDOWS=WINDOWS
											                        The size of smoothing windows if your inputfile is from scRNA-seq.
											                        The value is the number of genes which will be merge. Default value is 30.


-Permutation permutation, --Permutation=Permutation
												Performing reconstruction of permuted tree (T) or not (F) to estimate background distribution. If not, only permuted copy number profile will be used to perform LSA. Default value is F due to time cost.

```

Input files
===========

DNA input files:

	Two kinds of input files are allowed in MEDALT:

	(1) Integer copy number profile from scDNA-seq

	(2) Inferred integer copy number profile from scRNA-seq


  *scDNA-seq input*

  	chr	pos	cell1  cell2 cell3 ......
  	1	977836	2  3 1 ......
  	1	1200863	3 3 1	......

  *scRNA-seq input*

    	cell1  cell2 cell3 ......
    gene1	2  3 1 ......
    gene2	3 3 1	......

		For scRNA-seq data, we calculate an integer copy number by multiplying the relative CNA value from inferCNV by 2 (diploid) and rounding the results off to closest integers.

Run MEDALT package
============

    Python scTree.py [-O <output path>] [-W <smoothing window size>] –P <MEDALT package path> –I <input file> -D <input file type>
    [...] contains optional parameters. The mandatory arguments are -P, -I and -D. The input file type is either "D" or "R".
By default, -W is 30, which defines the smoothing window as 30 adjacent gene for scRNA-seq data.


Examples
========
Try MEDALT in the package directory on the different example datasets

**Example 1: Input integer copy number profile from scDNA-seq data**

	python scTree.py -P ./ -I ./example/scDNA.CNV.txt -D D -O ./example/output1

**Example 2: Input integer copy number profile from scRNA-seq data**

	python scTree.py -P ./ -I ./example/scRNA.inferCNV.txt -D R -O ./example/output2

>In order to save time, we don't reconstruct permuted trees. You can add it by -Permutation T

Output files
============

Three text files:

	(1) CNV.tree.txt file: an rooted directed tree includes three columns: parent node, child node and distance

	(2) segmental.LSA.txt file: significant broad fitness-associated CNAs

	(3) gene.LSA.txt file: significant focal fitness-associated CNAs.

Two figures:

	(1) singlecell.tree.pdf: the MEDALT of cells visualized by igraph. You can input CNV.tree.txt into Cytoscape to generate preferred visualization.

	(2) LSA.tree.pdf: the tree visualization of fitness-associated CNAs by igraph


Developer
=========
Fang Wang (fwang9@mdanderson.org), Qihan Wang (Chuck.Wang@rice.du)

Draft date
==========
April. 06, 2020
