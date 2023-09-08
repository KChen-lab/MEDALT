#main R function to perform linegae speciation analysis
#input: 1. the path of R script LSAfunction.R
#       2. input copy number profile from scDNA-seq or scRNA-seq
#       3. the infer tree file based on integer copy number profile
#       4. integer copy number matrix
#       5. output path of results
#       6. nuc_type "D" (DNA-seq) or "R" (RNA-seq)
#       7. genome version: hg19 or hg38
#       8. optional. If you have the permutation tree, input the path

suppressWarnings(suppressMessages({
    library(HelloRanges)
    library(igraph)
    library(stringr)
}))

args<-commandArgs(T)
PKG_PATH <- args[1]  # path to the package
CNV_PATH <- args[2]  # path to original raw input CNV
OUT_PATH <- args[3]  # path to output dir
SCT_PATH <- args[4]  # path to single cell tree
CNV_FILE <- args[5]  # path to dataTransfer-ed CNV
nuc_type <- args[6]  # nucleic acid type
ref__seq <- args[7]  # human genome version
if (length(args)==8){ permPATH=args[8] }
source(paste(PKG_PATH,"/SP2_LSA_FUNC.R",sep=""))
set.seed(1234) ; setwd(OUT_PATH)

REF_PATH <- paste0(PKG_PATH, "genomes/gencode_v", gsub("[^[:digit:]]", "", ref__seq), "_gene_pos.txt")
ONC_PATH <- paste0(PKG_PATH, "genomes/pathwaygene.txt")
BANDPATH <- paste(PKG_PATH, "genomes/", ref__seq, ".band.bed", sep="")


### single cell tree retrieval ###############################################################
### retrieve single cell tree inferred from the SC1_py_sctree.py                           ###
### then visualize the graph and find nodes/cells of interest to do LSA                    ###
##############################################################################################
cat("\n#####################################################\n")
cat(  "### now running SC4_RR_LSA_tree.R                 ###")
cat("\n#####################################################\n\n")
cat("Visualizing MEDALT\n")
# read in edges
sct_edges=read.csv(SCT_PATH, sep="\t", colClasses=c("character", "character", "integer")) # read edges
sct_nodes=data.frame(id=union(sct_edges$stt, sct_edges$end), size=3, color="lightgreen")  # initialize nodes
rownames(sct_nodes) <- sct_nodes$id                                                       # set rownames to be cell_id 

# find the root
root_node <- setdiff(sct_edges$stt, sct_edges$end)
root_coor <- which(sct_nodes$id==root_node)
sct_nodes$color[root_coor]="black"
# make graph obj
sct_graph <- graph_from_data_frame(d=sct_edges, vertices=sct_nodes, directed=T) # igraph obj

# use nodes' distance from root to color them
sct_nodes$dist <- as.integer(t(distances(sct_graph, v=V(sct_graph)[root_coor], to=V(sct_graph), weights=sct_edges$len)))
sct_nodes$dpth <- as.integer(t(distances(sct_graph, v=V(sct_graph)[root_coor], to=V(sct_graph), weights=rep(1, nrow(sct_edges)))))
rt_pal <- colorRampPalette(c("#fb8072", "#80b1d3"))
rt_col <- c("black", rt_pal(max(sct_nodes$dist))) ; sct_nodes$dist_col <- rt_col[sct_nodes$dist+1]
rt_col <- c("black", rt_pal(max(sct_nodes$dpth))) ; sct_nodes$dpth_col <- rt_col[sct_nodes$dpth+1]

# use edge "betweenness" to set the width of edges
eg_btw <- edge.betweenness(sct_graph)/20
eg_pal <- colorRampPalette(c("#3B3B3B", "#9C9C9C"))
eg_col <- eg_pal(max(sct_edges$len)+1)

graph_size <- 14
options(repr.plot.width=graph_size, repr.plot.height=graph_size)
sctlayout <- layout_with_kk(sct_graph, weights=(sct_edges$len+500))
pdf(file=paste(OUT_PATH, "/", "3_scTree_visualization.pdf", sep=""), width=graph_size, height=graph_size)  # , res=300, units="in", compression="lzw"
plot(sct_graph, layout=sctlayout,
     edge.arrow.size=.2, edge.width=sct_edges$dist+8, edge.label=sct_edges$len, edge.color="#7f7f7f",
     vertex.label=NA, vertex.color=sct_nodes$dpth_col, vertex.size=sct_nodes$size, 
     vertex.frame.width=3, vertex.frame.color="black")
try(dev.off(), silent=TRUE)

# calculate the depth, the number of children nodes for each cell
cell_pem <- data.frame(cells=sct_nodes$id, 
                       depth=sct_nodes$dpth,
                       child=sapply(sct_nodes$id, child, cellTree=sct_edges))

# perform lineage speciation analysis for cells that the number of children is no less than 5.
cell_pem <- cell_pem[cell_pem$child >= 5 & cell_pem$cells != root_node,] ### cells for permutation


### onco-gene & chr regions retrieval ########################################################
### from hg19/hg38, get the overlap of data CNV with onco-genes/chr-regions                ###
##############################################################################################
ref_gene <- read.csv(REF_PATH, sep="\t", header=F, colClasses=c("character", "character", "integer", "integer"))
colnames(ref_gene) <- c("gene", "chr", "stt", "end")
ref_gene <- ref_gene[grep("chr[^YM]", ref_gene$chr),]
#ref_gene$chr <- gsub("chr[xX]", "chr23", ref_gene$chr)

CNV_gene <- read.csv(CNV_PATH, sep="\t", header = TRUE)
if(!is.numeric(CNV_gene[,1])){
      rownames(CNV_gene)=CNV_gene[,1]
      CNV_gene = CNV_gene[,-1]}
CNV_gene <- CNV_gene[,intersect(sct_nodes$id, colnames(CNV_gene))]

cat("LSA segmentation!\n")
if (nuc_type=="D"){
    gene_loc     <- CNV_gene[,1:2]
    gene_loc[,1] <- paste("chr",gene_loc[,1],sep="")
    gene_loc$end <- gene_loc[,2]
    gene_loc <- gene_loc[!duplicated(gene_loc),]
    colnames(gene_loc) <- c("chr","stt","end")
}else if (nuc_type=="R"){
    CNV_gene <- round(CNV_gene)
    gene_idx <- match(row.names(CNV_gene), ref_gene$gene)    # subset to genes matched in ref_gene
    CNV_subs <- cbind(ref_gene[gene_idx[!is.na(gene_idx)], c("chr", "stt")], 
                      CNV_gene[         !is.na(gene_idx) ,]) # concat 2 dfs
    rownames(CNV_subs) <- rownames(CNV_gene)[!is.na(gene_idx)]
    
    gene_loc <- CNV_subs[,c(1,2,2)]
    colnames(gene_loc)=c("chr","stt","end")
    CNV_gene <- CNV_subs ; rm(CNV_subs)} ##################### CNV_gene is CNV_subs from now on

write.table(gene_loc, paste0(OUT_PATH, "/4_gene_loc.bed"), col.names=F, row.names=F, sep="\t", quote=F)  # for next bed-tool intersect

onc_pth <- read.csv(ONC_PATH, sep="\t")
colnames(onc_pth) <- c('chr', 'stt', 'end', 'name', 'pathway')
onc_pth <- onc_pth[onc_pth$name %in% ref_gene$gene,]             # onco genes in data
idx_onc <- match(onc_pth$name, ref_gene$gene)                    # correct gene location
onc_pth[,c("chr", "stt", "end")] <- ref_gene[idx_onc,c("chr", "stt", "end")]

#Gene level copy number profile
if (nuc_type=="D"){
    CNV__DNA <- read.csv(CNV_FILE,sep="\t")
    CNV_onco <- do.call(cbind,lapply(1:dim(onc_pth)[1], geneCNAfunction, onc_pth, CNV__DNA, gene_loc))
    colnames(CNV_onco) <- as.character(onc_pth$name)
    rownames(CNV_onco) <- rownames(CNV__DNA)
}else if (nuc_type == "R"){
    CNV_gene <- CNV_gene[, 3:dim(CNV_gene)[2]]  # get rid of chr & loc, remaining only cell barcodes in columns
    idx_onco <- match(onc_pth$name, rownames(CNV_gene))
    CNV_onco <- t(CNV_gene[idx_onco[!is.na(idx_onco)],])
    CNV_onco <- round(CNV_onco)
    colnames(CNV_onco) <- onc_pth$name[!is.na(idx_onco)]}

index_NA <- apply(CNV_onco, 2, find_NAs)    # get rid off NAs
CNV_onco <- CNV_onco[, index_NA==1]
rownames(CNV_onco) <- str_replace_all(rownames(CNV_onco), "\\.", "_")

ref_band = read.csv(BANDPATH, sep="\t", header = F, 
                    colClasses=c("character", "integer", "integer", "character"))
colnames(ref_band) <- c("chr", "stt", "end", "band")

chrom = unique(ref_band[,1])
chr__arm = c()
band_reg = c()
for (i in chrom){                                    # get band information of all chr
    sub__ref <- ref_band[ref_band[,1]==i,]   
    chr_parm <- sub__ref[grep("p",sub__ref[,4]),]    # find the p arm region
    chr_qarm <- sub__ref[grep("q",sub__ref[,4]),]    # find the q arm region
    arm_coor <- data.frame(chr=i,                    # put p/q arm info to a dataframe
                         start=c(head(chr_parm,1)[,2], head(chr_qarm,1)[,2]), 
                           end=c(tail(chr_parm,1)[,3], tail(chr_qarm,1)[,3]), band=c("p","q"))
    chr__arm <- rbind(chr__arm, arm_coor)
    sub.band <- do.call(rbind, strsplit(sub__ref[,4], split="[.]"))
    bandname <- unique(sub.band[,1])
    arm_posi <- sapply(bandname, get_arm_position, sub.band, sub__ref)
    sub_regn <- data.frame(chr=i, start=arm_posi[1,], end=arm_posi[2,], band=bandname)
    band_reg <- rbind(band_reg, sub_regn)}

chr__arm$length <- chr__arm$end - chr__arm$start
band_reg$length <- band_reg$end - band_reg$start
ref_band$ID     <- paste(ref_band$chr, ref_band$band, sep=":")

bed_command = bedtools_intersect(paste0("-a ", PKG_PATH, "/genomes/", ref__seq, ".band.bed ", 
                                        "-b ", paste0(OUT_PATH, "/4_gene_loc.bed")))
bed_outs <- as.data.frame( eval(bed_command) )
bed_outs$ID  <- paste(bed_outs$seqnames, ":", bed_outs$name, sep="") # band id: chr_band
bed_outs$ID1 <- paste(bed_outs$seqnames, "_", bed_outs$end,  sep="") # band id: chr_location
gene_loc$ID1 <- paste(gene_loc$chr, "_", gene_loc$stt, sep="")
bed_outs <- bed_outs[!duplicated(bed_outs),]
ID=unique(bed_outs$ID)

### get CN matrix for chr band
CNV_band = do.call(cbind, lapply(ID, subsetCNV, CNV_gene=CNV_gene, bed_outs=bed_outs, gene_loc=gene_loc))
colnames(CNV_band)=ID

### perform LSA ##############################################################################
##############################################################################################
### perform LSA analysis
geneGscore <- lapply(cell_pem$cells, lineageScore, CNV_onco, sct_edges)
names(geneGscore) <- cell_pem$cells

#calculte CFL at genomic bin level
bandGscore=lapply(cell_pem$cells, lineageScore, CNV_band, sct_edges)
names(bandGscore)=cell_pem$cells

### run permutation to estimate significance of LSA results ##################################
##############################################################################################
# do calculation for permutation dataset
cat("Calculating permutation CFL\n")
t0 <- Sys.time() ; tt <- Sys.time()
arg_8 <- length(args) < 8
if (arg_8){ # if no permutation tree, calculate CFL based on real tree structure
    times=200 ; rpfeq=times%/%10
    permuteres=lapply(1:times, denovo_perm_score, 
                      CNV_gene, ID, bed_outs, nuc_type, onc_pth, gene_loc, ref_gene, sct_edges)
} else {    # permutation trees exist, calculate CFL based on permuted tree structure
    permutefile <- list.files(permPATH)
    permutefile <- permutefile[grep("celltree", permutefile)]
    cat(paste("There are ", length(permutefile), " permutation trees."))
    if (length(permutefile)>0){
        permuteres = lapply(permutefile, permuteTreeScore, ID, bed_outs, nuc_type, onc_pth, gene_loc, ref_gene, permPATH)}}

cat(paste("runtime:", Sys.time()-t0, "\n"))

#Estimate emperical p value
cat("Estimate emperical p value\n")
p_value <- lapply(1:nrow(cell_pem), significanceLevel, bandGscore, geneGscore, permuteres, cell_pem)
pcutoff <- 0.01
gene_sig <- CollectAsso(p_value, pcutoff, sct_edges, cell_pem)$generes # at genomic bin and gene level
band_sig <- CollectAsso(p_value, pcutoff, sct_edges, cell_pem)$bandres # estimate emperical p value
gene_sig <- unique(gene_sig)
band_sig <- do.call(rbind, lapply(unique(as.character(band_sig$cell)), mergeCNA, band_sig, band_reg, chr__arm, ref_band))
band_sig <- unique(band_sig)

LSAres=list()
if (!is.null(gene_sig)){
    index <- match(gene_sig$cell, cell_pem$cells)
    gene_sig$child <- cell_pem$child[index]
    LSAres$geneLSA <- gene_sig
    paraGene <- table(as.character(gene_sig$region))
    paraGene <- paraGene[paraGene>1] #gene associated with more than one independent lieage
    if (length(paraGene)>0){
        paraGenesig <- GenePara(gene_sig, permuteres, type="gene", cell_pem)#parallele evolution estimation
        if (!is.null(paraGenesig)){ LSAres$paraGene <- paraGenesig }}}

if (!is.null(band_sig)){
    index <- match(band_sig$cell, cell_pem$cells)
    band_sig$child <- cell_pem$child[index]
    LSAres$bandLSA <- do.call(rbind, lapply(unique(band_sig$cell), CombineRegion, band_sig, ref_band))
    paraBand <- table(as.character(band_sig$region))
    paraBand <- paraBand[paraBand>1] #CNAs of genomic bin associated with >1 independent lineages
    if (length(paraBand)>0){
        paraBandsig <- GenePara(band_sig, permuteres, type="band", cell_pem) #parallel evolution test
        if (!is.null(paraBandsig)){ LSAres$paraBand <- paraBandsig }}}

cat("Estimate parallel evolution\n")
all_LSA=c()
if ("geneLSA" %in% names(LSAres)){
    geneLSA <- LSAres$geneLSA
    geneLSA$CNA[geneLSA$Score>0] <- "AMP"
    geneLSA$CNA[geneLSA$Score<0] <- "DEL"
    all_LSA <- rbind(all_LSA, geneLSA)
    write.table(geneLSA, paste(OUT_PATH, "/5_LSA_gene.tsv", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}else{ cat("No LSA is identified at gene level!\n") }

if ("bandLSA" %in% names(LSAres)){
    bandLSA=LSAres$bandLSA
    bandLSA$CNA[bandLSA$Score>0]="AMP"
    bandLSA$CNA[bandLSA$Score<0]="DEL"
    all_LSA <- rbind(all_LSA, bandLSA)
    write.table(bandLSA, paste(OUT_PATH,"/5_LSA_band.tsv", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}else{ cat("No segmental LSA is identified!\n") }

paraEvent=c()
if ("paraBand" %in% names(LSAres)){ paraEvent=rbind(paraEvent,LSAres$paraBand) }
if ("paraGene" %in% names(LSAres)){ paraEvent=rbind(paraEvent,LSAres$paraGene) }
if (!is.null(paraEvent)){
  paraEvent <- paraEvent[!is.na(paraEvent$pvalue),]
  if (dim(paraEvent)[1]!=0){
    write.table(paraEvent, paste(OUT_PATH,"/5_LSA_parallel.tsv",sep=""), col.names=T, row.names=F, sep="\t", quote=FALSE) }}

### embed LSA result back to scTree for result interpretation ################################
##############################################################################################
LSA_to_plot = list()
if ("geneLSA" %in% ls()) {LSA_to_plot[["geneLSA"]] <- geneLSA}
if ("bandLSA" %in% ls()) {LSA_to_plot[["bandLSA"]] <- bandLSA}
if ("all_LSA" %in% ls()) {LSA_to_plot[["all_LSA"]] <- all_LSA}

for (LSA_name in names(LSA_to_plot)){
    LSA_use <- LSA_to_plot[[LSA_name]]
    sct_nodes$CNA_gene_count <- 0         # set 0 to all CNA_gene_counts
    lsa_n_cna <- table(LSA_use$cell)      # n_cna in each lsa node
    sct_nodes[names(lsa_n_cna), "CNA_gene_count"] <- lsa_n_cna

    sct_nodes$CNA_gene_annot <- ""
    for (cell_i in names(lsa_n_cna)){
        CNA <- LSA_use[LSA_use$cell==cell_i,]                               # get CNAs for single cells
        CNA <- CNA[order(-sign(CNA$Score)*log(CNA$adjustp)),]               # sort CNAs by significance & change direction
        nl_coor <- seq(1,nrow(CNA),3)                                       # CNA counts by 3
        CNA$CNA[nl_coor] <- gsub("^", "\n", CNA$CNA[nl_coor])               # insert new lines
        cna <- paste(CNA$CNA, str_pad(CNA$region, 7, "left", "_"), sep="_", collapse=" ; ")   # concatenate all CNAs
        sct_nodes[cell_i, "CNA_gene_annot"] <- paste0(cell_i, ":", nrow(CNA), cna)}  # record to single cell tree graph
    sct_nodes[root_node, "CNA_gene_count"] <- max(sct_nodes$CNA_gene_count)
    
    cna_edges <- sct_edges[sct_edges$end %in% cell_pem$cells,]
    cna_nodes <- sct_nodes[union(cna_edges$stt, cna_edges$end),]
    cna_graph <- graph_from_data_frame(d=cna_edges, vertices=cna_nodes, directed=T)

    # plot CNA graph
    plot_size <- 16
    try(dev.off(), silent=TRUE)
    pdf(file=paste(OUT_PATH, "/5_scTree_", LSA_name, "_CNA_RT.pdf", sep=""), width=graph_size, height=graph_size)  # , res=300, units="in", compression="lzw"
    plot(cna_graph, layout=layout_as_tree, 
        edge.arrow.size=.2, edge.width=cna_edges$dist+8, edge.label=cna_edges$dist, edge.color="#9c9c9c",
        vertex.frame.width=3, vertex.frame.color="black", vertex.label=cna_nodes$CNA_gene_annot, vertex.color=cna_nodes$color, vertex.size=sqrt(cna_nodes$CNA_gene_count)*2)
    try(dev.off(), silent=TRUE)

    # plot CNA graph
    plot_size <- 16
    options(repr.plot.width=plot_size, repr.plot.height=plot_size)
    pdf(file=paste(OUT_PATH, "/5_scTree_", LSA_name, "_CNA_KK_cna.pdf", sep=""), width=graph_size, height=graph_size)  # , res=300, units="in", compression="lzw"
    plot(cna_graph, layout=sctlayout[match(union(cna_edges$stt, cna_edges$end), rownames(sct_nodes)),], 
        edge.arrow.size=.2, edge.width=cna_edges$dist+8, edge.label=cna_edges$dist, edge.color="#9c9c9c",
        vertex.frame.width=3, vertex.frame.color="black", vertex.label=cna_nodes$CNA_gene_annot, vertex.color=cna_nodes$color, vertex.size=sqrt(cna_nodes$CNA_gene_count)*2)
    try(dev.off(), silent=TRUE)
    
    # plot CNA graph
    plot_size <- 16
    options(repr.plot.width=plot_size, repr.plot.height=plot_size)
    pdf(file=paste(OUT_PATH, "/5_scTree_", LSA_name, "_CNA_KK_sct.pdf", sep=""), width=graph_size, height=graph_size)  # , res=300, units="in", compression="lzw"
    plot(sct_graph, layout=sctlayout, 
        edge.arrow.size=.2, edge.width=sct_edges$dist+8, edge.label=sct_edges$dist, edge.color="#9c9c9c",
        vertex.frame.width=3, vertex.frame.color="black", vertex.label=sct_nodes$CNA_gene_annot, vertex.color=sct_nodes$color, vertex.size=sqrt(sct_nodes$CNA_gene_count)*2)
    try(dev.off(), silent=TRUE)




}