#R function for permutation process
#input: 1. the path of R script LSAfunction.R
#       2. copy number profile matrix
#       3. datatype "D" (DNA-seq) or "R" (RNA-seq)
#       4. set the path which store the tree corresponding to permutation datasets
#       5. the number of genes to estimate copy number of genomic bin if you set datatype as R
#           defaul 30 genes

args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
datatype=args[3]
permutationPath=args[4]
if (length(args)==6){
  delt=args[5]
  hg=args[6]
  if (hg == "hg19"){
    reference=read.csv(paste(datapath,"/gencode_v19_gene_pos.txt",sep=""),sep="\t",header=F)
  }else{
    reference=read.csv(paste(datapath,"/gencode_v38_gene_pos.txt",sep=""),sep="\t",header=F)

  }
}
set.seed(1234)
source(paste(datapath,"/LSAfunction.R",sep=""))
data = read.csv(inputfile,sep="\t",header = TRUE)

#refernce genome to check if input data are ordered
#this is only used for RNA-seq data
if (datatype=="R"){
  data=round(data*2)#integer copy number
}
for (j in 1:100){

  #permute copy number profile from DNA-seq data
  if (datatype=="D"){
    region=data[,1:2]
    region[,1]=paste("chr",region[,1],sep="")
    region$end=region[,2]
    colnames(region)=c("chrom","chrompos","end")
    region$ID1=paste(region$chrom,"_",region$chrompos,sep="")
    permuteCNV=permuteSeg(data,region)
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }else if (datatype=="R"){

    #Permutate the copy number profile by genes from the same chromosome into different cells
    permuteCNV=permuteGene(data,reference)

    #calculate copy number of genomic bin
    regionCNV=BinCNV(reference,permuteCNV,as.numeric(delt))
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".gene.CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
    write.table(t(regionCNV),paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }
}
