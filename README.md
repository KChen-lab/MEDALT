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
  -R PERMUTATION, --PERMUTATION=PERMUTATION
                        Performing tree reconstruction based on permutation data (T) or not (F) to estimate background distribution.
                        If yes, both permuted tree and copy number profile will be used to perform LSA. Otherwise, only permuted copy number profile will be used to perform LSA.
                        Default value is F due to time cost.

```

Input files
===========

Single cell copy number input files:

	Two kinds of input files are allowed in MEDALT:

	(1) Integer copy number profile from scDNA-seq

	(2) Inferred copy number profile from scRNA-seq

  *scDNA-seq input*

  	chr	pos	cell1  cell2 cell3 ......
  	1	977836	2  3 1 ......
  	1	1200863	3 3 1	......

  *scRNA-seq input*

    	cell1  cell2 cell3 ......
    gene1	0.5  1.5 2.1 ......
    gene2	1.1 1.8 0.6	......

>For scRNA-seq input, the copy number is inferred relative copy number (relative to normal cells) instead of integer copy number. If value = 1, it means diploid. We directly incorporate inferCNV result (expr.infercnv.12_denoised.dat) as input.

Run MEDALT package
============

    Python scTree.py [-O <output path>] [-W <smoothing window size>] [-R <permutation tree reconstruction>] –P <MEDALT package path> –I <input file> -D <input file type>
    [...] contains optional parameters.
		The mandatory arguments are -P, -I and -D.
		The input file type (-D) is either "D" (DNA) or "R" (RNA).
>By default, we estimate background using by-chromosome permuted single cell copy number matrix rather than reconstructing a tree from permuted matrix due to time cost. You can change the setting by -R T. The default value of smoothing window size (-W) is 30, which defines the smoothing window as 30 adjacent genes for scRNA-seq data.  


Examples
========
Try MEDALT in the package directory on the different example datasets

**Example 1: Input integer copy number profile from scDNA-seq data**

	python scTree.py -P ./ -I ./example/scDNA.CNV.txt -D D -O ./example/outputDNA

**Example 2: Input inferred relative copy number profile from scRNA-seq data**

	python scTree.py -P ./ -I ./example/scRNA.inferCNV.txt -D R -O ./example/outputRNA

>In order to save time, we don't reconstruct trees based on permutation data. You can set -R T
to reconstruct permuted tree.

Output files
============

Three text files:

	(1) CNV.tree.txt which is an rooted directed tree including three columns: parent node, child node and distance

	(2) segmental.LSA.txt which includes broad CNAs significantly associated with lineage expansion

	(3) gene.LSA.txt which includes focal (gene) CNAs significantly associated with lineage expansion

> If there is parallel evolution event, the results will be saved in a separate file.

Two figures:

	(1) singlecell.tree.pdf which is a visualization of MEDALT by igraph. You also can input CNV.tree.txt into Cytoscape to generate preferred visualization.

	(2) LSA.tree.pdf which is a visualization of identified CNAs by igraph.

> In LSA figure, we only show top 3 events for each lineage. You can check more details in segmental or gene level LSA file.


Developer
=========
Fang Wang (fwang9@mdanderson.org), Qihan Wang (Chuck.Wang@rice.du)

Draft date
==========
April. 06, 2020
