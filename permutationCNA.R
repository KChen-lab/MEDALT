#R function for permutation process
#input: 1.the path of R script LSAfunction.R
#       2.copy number profile matrix
#       3.datatype "D" (DNA-seq) or "R" (RNA-seq)
#       4.set the path which store the tree corresponding to permutation datasets
#       5.the number of genes 
args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
datatype=args[3]
permutationPath=args[4]
if (length(args)==5){
  delt=args[5]
}
set.seed(1234)
source(paste(datapath,"/LSAfunction.R",sep=""))
data = read.csv(inputfile,sep="\t",header = TRUE)
reference=read.csv(paste(datapath,"/gencode_v19_gene_pos.txt",sep=""),sep="\t",header=F)
if (datatype=="R"){
  data=round(data*2)
}
for (j in 1:100){
  if (datatype=="D"){
    region=data[,1:2]
    region[,1]=paste("chr",region[,1],sep="")
    region$end=region[,2]
    colnames(region)=c("chrom","chrompos","end")
    region$ID1=paste(region$chrom,"_",region$chrompos,sep="")
    permuteCNV=permuteSeg(data,region)
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }else if (datatype=="R"){
    permuteCNV=permuteGene(data,reference)
    regionCNV=BinCNV(reference,permuteCNV,as.numeric(delt))
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".gene.CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
    write.table(t(regionCNV),paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }
}
