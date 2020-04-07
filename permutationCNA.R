args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
datatype=args[3]
permutationPath=args[4]
delt=args[5]
source(paste(datapath,"/LSAfunction.R",sep=""))
data = read.csv(inputfile,sep="\t",header = TRUE)
reference=paste(datapath,"/gencode_v19_gene_pos.txt",sep="")
for (j in 1:100){
  if (datatype=="D"){
    permuteCNV=permuteSeg(data)
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }else if (datatype=="R"){
    permuteCNV=permuteGene(data,reference)
    regionCNV=BinCNV(refgene,permuteCNV,delt)
    write.table(permuteCNV,paste(permutationPath,"/permute.",j,".gene.CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
    write.table(regionCNV,paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),col.names = T, row.names=T,quote=F,sep="\t")
  }
}
