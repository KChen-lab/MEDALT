####transfer format of input data
##Input integer copy number profile with header from scDNA-seq, separeted by tab.
##The row in DNA input data is chromosomal bin. The range of chromosome is from 1 to 23.
##The first column is chromosome, the second column is the start position of bin, and the following columns are sequenced cells.
DNAinput <- function(inputfile){
  data=read.csv(inputfile,sep="\t")
  region=table(data[,1])
  regionseq=c()
  for (i in 1:length(region)){
    regionseq=c(regionseq,c(1:region[i]))
  }
  row.names(data)=paste("chr",data$chrom,"_",regionseq,sep="")
  CNV=data[,3:dim(data)[2]]
  CNV=t(CNV)
  write.table(CNV, paste(inputfile,".CNV.txt",sep=""), sep = "\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  #return(paste(inputfile,".CNV.txt",sep=""))
}

####transfer format of input data
##Input the expression data inferred from InferCNV: infercnv_obj@expr.data
##input the reference for gene position information: geneName chromosome  start end
##delt means the size of bin which is defined by the number of gene
RNAinput <- function(inputfile, reference, delt){
  data=read.csv(inputfile,sep="\t")
  geneInfo=read.csv(reference,sep="\t",header=F)
  index=match(row.names(data),as.character(geneInfo[,1]))
  newdata=cbind(geneInfo[index[!is.na(index)],2:4],data[!is.na(index),])
  ll=nchar(as.character(newdata[,1]))
  chromo=data.frame(chr=as.character(newdata[,1]),ll=ll)
  chro=apply(chromo,1,function(x){
    return(substr(x[1],start=4,stop=x[2]))
  })
  chro[chro=="X"]=23
  newdata[,1]=chro
  newdata=newdata[newdata[,1]!="M"&newdata[,1]!="Y",]
  chrom=c(1:23)
  chrom=intersect(chrom,chro)
  segdata=c()
  chrregion=c()
  for (i in chrom){
    subseg=c()
    subdata=newdata[newdata[,1]==i,]
    subdata=subdata[order(as.numeric(as.character(subdata[,2]))),]
    kk=dim(subdata)[1]/delt
    intekk=round(kk)
    if (intekk >1){
      for (j in 1:(intekk-1)){
        sub1=subdata[((j-1)*delt+1):(j*delt),4:dim(subdata)[2]]
        subseg=rbind(subseg,apply(sub1,2,mean))
        chrregion=c(chrregion,paste(i,"_",j,sep=""))
      }
      subseg=rbind(subseg,apply(subdata[((intekk-1)*delt+1):dim(subdata)[1],4:dim(subdata)[2]],2,mean))
      chrregion=c(chrregion,paste(i,"_",intekk,sep=""))
    }else{
      subseg=apply(subdata[,4:dim(subdata)[2]],2,mean)
      chrregion=c(chrregion,paste(i,"_",1,sep=""))
    }
    segdata=rbind(segdata,subseg)
  }
  row.names(segdata)=paste("chr",chrregion,sep="")
  segdata=t(round(segdata*2))
  write.table(segdata, paste(inputfile,".CNV.txt",sep=""), sep = "\t",col.names = TRUE, row.names = TRUE, quote = FALSE)
  #return(paste(inputfile,".CNV.txt",sep=""))
}

args<-commandArgs(T)
inputfile=args[1]
datatype=args[2]
if (datatype=="D"){
  DNAinput(inputfile)
}else if (datatype == "R"){
  datapath=args[3]
  delt = args[4]
  RNAinput(inputfile,reference=paste(datapath,"/gencode_v19_gene_pos.txt",sep=""),delt=as.numeric(delt))
}
