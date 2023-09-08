####transfer format of input data
##Input integer copy number profile with header from scDNA-seq, separeted by tab.
##The row in DNA input data is chromosomal segment. The range of chromosome is from 1 to 23.
##The first column is chromosome, the second column is the start position of bin, and the following columns are sequenced cells.
DNAinput <- function(in_file, outfile){
    cnv_mtx <- read.csv(in_file, sep="")
    regions <- table(cnv_mtx[,1])
    rgn_seq <- c()
    for (i in 1:length(regions)){
        rgn_seq=c(rgn_seq,c(1:regions[i]))
    }
    row.names(cnv_mtx)=paste("chr",cnv_mtx[,1],"_",rgn_seq,sep="")
    CNV=cnv_mtx[,3:dim(cnv_mtx)[2]]
    CNV=t(CNV)
    write.table(CNV, outfile, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    #return(paste(in_file,".CNV.txt",sep=""))
}

#### transfer format of input data
## Input the expression data inferred from InferCNV: infercnv_obj@expr.data
## input the reference for gene position information: geneName chromosome start end
## bin_siz means the size of bin which is defined by the number of gene
RNAinput <- function(in_file, ref_pth, bin_siz, outfile){
    cnv_mtx=read.csv(in_file,sep="") # "" to cope with any kinds of \t, \s
    if(!is.numeric(cnv_mtx[,1])){      # if the first col is gene name
        rownames(cnv_mtx)=cnv_mtx[,1]  # make it row names
        cnv_mtx = cnv_mtx[,-1]}        # then remove it from the mtx
    
    cnv_mtx <- round(cnv_mtx)          # round cn profile
    ref_gen <- read.csv(ref_pth, sep="\t", header=F)                # reference genes
    idx_gen <- match(row.names(cnv_mtx), as.character(ref_gen[,1])) # use only matched genes
    newdata <- cbind(ref_gen[idx_gen[!is.na(idx_gen)],2:3],         # chr & start
                     cnv_mtx[        !is.na(idx_gen) ,])
    cat(paste0("", dim(newdata)[1], "/", dim(cnv_mtx)[1], " genes matched in ref_seq.\n"))
    
    newdata[,1] <- gsub("chr", "", as.character(newdata[,1]))
    #newdata[,1] <- gsub("X", "23", newdata[,1])
    newdata <- newdata[newdata[,1]!="M" & newdata[,1]!="Y",]
    
    seg_dat <- c() ; chr_rgn <- c()
    for (i in unique(newdata[,1])){                          # within each chr
        sub_seg <- c()
        subdata <- newdata[newdata[,1]==i,]
        subdata <- subdata[order(as.numeric(subdata[,2])),]  # sort by chr coor
        bin_num <- round(dim(subdata)[1]/bin_siz)            # n bins in this chr
        if(bin_num>0){
            # take the mean for 1:(n-1) bin
            for (j in 1:max(1, bin_num-1)){
                sub_bin <- subdata[((j-1)*bin_siz+1):(j*bin_siz), -1:-2] 
                sub_seg <- rbind(sub_seg, apply(sub_bin, 2, mean))
                chr_rgn <- c(chr_rgn, paste(i,"_", j, sep=""))}
            # take the mean for last bin
            lastbin <- ((bin_num-1)*bin_siz+1):dim(subdata)[1] # mean last bin
            sub_seg <- rbind(sub_seg, apply(subdata[lastbin,-1:-2],2,mean))
            chr_rgn <- c(chr_rgn, paste(i, "_", bin_num, sep=""))
            # attach to overall data
            seg_dat <- rbind(seg_dat, sub_seg)
        }else{
            sub_seg <- apply(subdata[,-1:-2],2,mean)
            chr_rgn <- c(chr_rgn, i)
            seg_dat <- rbind(seg_dat, sub_seg)
    }}
    
    row.names(seg_dat) <- paste("chr", chr_rgn, sep="")
    seg_dat=t(round(seg_dat))  # eventually wants to make it integer
    write.table(seg_dat, outfile, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    cat(paste0("saved file: ", outfile, "\n"))}

#################################################################
## main run                                                    ##
#################################################################
cat("\n#####################################################\n")
cat(  "### now running SC2_RR_dataTransfer.R             ###")
cat("\n#####################################################\n\n")
cmd_arg <- commandArgs(T)
outpath <- cmd_arg[1]
in_file <- cmd_arg[2]
nt_type <- cmd_arg[3]
outfile <- cmd_arg[4]
ref_pth <- cmd_arg[5]
if (nt_type == "R"){ bin_siz <- cmd_arg[6] }
setwd(outpath)

if (nt_type=="D"){
    DNAinput(in_file, outfile)
}else if (nt_type == "R"){
    bin_siz <- as.numeric(bin_siz)
    RNAinput(in_file, ref_pth=ref_pth, bin_siz=bin_siz, outfile=outfile)
}else{
    cat("invalid specified datatype, should be either D or R\n")}