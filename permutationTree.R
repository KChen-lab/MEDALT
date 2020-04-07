args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
datatype=args[3]
permutationPath=args[4]

data = read.csv(inputfile,sep="\t",header = TRUE)
reference=paste(datapath,"/gencode_v19_gene_pos.txt",sep="")
if (datatype=="D"){
  region=data[,1:2]
  region[,1]=paste("chr",region[,1],sep="")
  region$end=region[,2]
  colnames(region)=c("chrom","chrompos","end")
}else if (datatype=="R"){
  geneInfo=read.csv(reference,sep="\t",header=F)
  index=match(row.names(data),as.character(geneInfo[,1]))
  newdata=cbind(geneInfo[index[!is.na(index)],2:3],data[!is.na(index),])
  rownames(newdata)=rownames(data)[!is.na(index)]
  newdata=newdata[as.character(newdata[,1])!="chrM"&as.character(newdata[,1])!="chrY",]
  newdata[,1]=as.character(newdata[,1])
  newdata[newdata[,1]=="chrX",1]="chr23"
  region=newdata[,1:2]
  region[,3]=region[,2]
  colnames(region)=c("chrom","chrompos","end")
  cnv=t(data)
  data=newdata
}
times=100
if (datatype == "D"){
  permuteCNV=permuteSeg(cnv)
}
if (datatype == "R"){
  permuteCNV=permuteGene(cnv,reference)
  
}
Gscore=lapply(as.character(cell1$cell),lineageScore,permuteCNV1,celltree)
names(Gscore)=as.character(cell1$cell)
geneGscore=lapply(as.character(cell1$cell),lineageScore,geneCNV,celltree)
names(geneGscore)=as.character(cell1$cell)
return(list(bandGscore=Gscore,geneGscore=geneGscore))
