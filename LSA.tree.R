r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
if (!"HelloRanges" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("HelloRanges")
}
if (!"igraph" %in% installed.packages()){
  install.packages("igraph")
}
if (!"DescTools" %in% installed.packages()){
  install.packages("DescTools")
}
library(HelloRanges)
library(igraph)
library(DescTools)
args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
treeName=args[3]
CNVfile=args[4]
outpath = args[5]
datatype=args[6]
if (length(args)==7){
  permutationPath=args[7]
}
source(paste(datapath,"/LSAfunction.R",sep=""))
######Plot cell tree
print.noquote("Visualization MEDALT!")
celltree=read.csv(treeName,sep="\t")
nodes=data.frame(id=union(as.character(celltree[,1]),as.character(celltree[,2])),size=5)
nodes$color="lightblue"
nodes$color[nodes$id==setdiff(as.character(celltree[,1]),as.character(celltree[,2]))]="black"
net <- graph_from_data_frame(d=celltree, vertices=nodes, directed=T)
pdf(file=paste(outpath,"/singlecell.tree.pdf",sep=""),width = 5,height = 5,useDingbats = F)
plot(net, vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label=NA)
dev.off()
######input
data = read.csv(inputfile,sep="\t",header = TRUE)
reference=paste(datapath,"/gencode_v19_gene_pos.txt",sep="")
refer.band=read.csv(paste(datapath,"/band.bed",sep=""),sep="\t",header = F)
chrom=as.character(unique(refer.band[,1]))
arm.band=c()
band.region=c()
for (i in chrom){
  subrefer=refer.band[refer.band[,1]==i,]
  parm=subrefer[grep("p",subrefer[,4]),]
  qarm=subrefer[grep("q",subrefer[,4]),]
  armregion=data.frame(chr=i,start=c(parm[1,2],qarm[1,2]),end=c(parm[dim(parm)[1],3],qarm[dim(qarm)[1],3]),band=c("p","q"))
  arm.band=rbind(arm.band,armregion)
  subband=do.call(rbind,strsplit(as.character(subrefer[,4]),split="[.]"))
  bandname=unique(subband[,1])
  pos=sapply(bandname,function(x,subband,subrefer){
    k=which(subband[,1]==x)
    start=subrefer[min(k),2]
    end=subrefer[max(k),3]
    return(c(start,end))
  },subband,subrefer)
  subregion=data.frame(chr=i,start=pos[1,],end=pos[2,],band=bandname)
  band.region=rbind(band.region,subregion)
}
arm.band$length=arm.band$end-arm.band$start
band.region$length=band.region$end-band.region$start
refer.band$ID=paste(as.character(refer.band$V1),as.character(refer.band$V4),sep=":")

###segmentation
print.noquote("LSA segmentation!")
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
write.table(region,"region.bed",col.names = F,row.names = F,sep="\t",quote = F)
code=bedtools_intersect(paste("-a ",datapath,"/band.bed -b region.bed",sep=""))
ans <- eval(code)
ans=as.data.frame(ans)
ans$ID=paste(ans$seqnames,":",ans$name,sep="")
ans$ID1=paste(ans$seqnames,"_",ans$end,sep = "")
region$ID1=paste(region$chrom,"_",region$chrompos,sep="")
ID=unique(ans$ID)
newCNV=do.call(cbind,lapply(ID, function(id,data,ans,region){
  chrID=ans$ID1[ans$ID==id]
  index=match(chrID,region$ID1)
  if (length(index)==1){
    return(as.numeric(data[index,3:dim(data)[2]]))
  }else{
    return(round(apply(data[index,3:dim(data)[2]],2,mean)))
  }
},data=data,ans=ans,region))
colnames(newCNV)=ID
####lineage partitioning to define CFL
print.noquote("Calculating CFL")
cell=union(as.character(celltree[,1]),as.character(celltree[,2]))
cell=data.frame(cell=cell)
cell$depth=sapply(as.character(cell$cell),depthFunction,cellTree=celltree)
cell$subtreesize=sapply(as.character(cell$cell),subtreeSize,cellTree=celltree)
cell1=cell[cell$subtreesize>=5,]
cell1=cell1[cell1$cell!="root",]
Gscore=lapply(as.character(cell1$cell),lineageScore,newCNV,celltree)
names(Gscore)=as.character(cell1$cell)
pathwaygene=read.csv(paste(datapath,"/pathwaygene.txt",sep=""),sep="\t")
if (datatype=="D"){
  cnv=read.csv(CNVfile,sep="\t")
  oncogenicCNV=do.call(cbind,lapply(1:dim(pathwaygene)[1],geneCNAfunction,pathwaygene=pathwaygene,ancestorCNV=cnv,generegion=region))
  colnames(oncogenicCNV)=as.character(pathwaygene$name)
}else if (datatype == "R"){
  data=data[,3:dim(data)[2]]
  index=match(as.character(pathwaygene$name),rownames(data))
  oncogenicCNV=t(data[index[!is.na(index)],])
  oncogenicCNV=round(oncogenicCNV)
  colnames(oncogenicCNV)=as.character(pathwaygene$name)[!is.na(index)]
}
index=apply(oncogenicCNV,2,function(x){
  if (NA %in% x){
    return(0)
  }else{
    return(1)
  }
})
oncogenicCNV=oncogenicCNV[,index==1]
geneGscore=lapply(as.character(cell1$cell),lineageScore,oncogenicCNV,celltree)
names(geneGscore)=as.character(cell1$cell)
realres=list(cell=cell1,bandGscore=Gscore,geneGscore=geneGscore)
#############permutation
print.noquote("Calculating permutation CFL")
if (length(args) < 7){
  times=100
  permuteres=lapply(1:times,function(j,data,ID,ans,datatype,pathwaygene,generegion,reference){
    score=permuteScore(data,ID,ans,datatype,pathwaygene,generegion=region,reference)
    return(score)
  },data,ID,ans,datatype,pathwaygene,generegion=region,reference)
}else if (length(args)==7){
  permutefile=list.files(permutationPath)
  permutefile=permutefile[grep("celltree",permutefile)]
  times=length(permutefile)
  print.noquote(paste("There are ",length(permutefile)," permutation trees."))
  if (length(permutefile)>0){
    permuteres=lapply(permutefile,permuteTreeScore,ID,ans,datatype,pathwaygene,generegion=region,reference,permutationPath)
  }
}
#####Estimate emperical p value
print.noquote("Estimate emperical p value")
realcell=realres$cell
#index=match("root",as.character(realcell$cell))
realband=realres$bandGscore
realgene=realres$geneGscore
#if (!is.na(index)){
#  realband[[index]] <- NULL
#  realgene[[index]] <- NULL
#  realcell=realcell[realcell$cell!="root",]
#}
pvalue=lapply(1:dim(realcell)[1],significanceLevel,realband,realgene,permuteres,realcell)
bandsig=CollectAsso(pvalue,cutoff=1/times,celltree,realcell)$bandres
genesig=CollectAsso(pvalue,cutoff=1/times,celltree,realcell)$generes
bandsig=do.call(rbind,lapply(unique(as.character(bandsig$cell)),mergeCNA,bandsig,band.region,arm.band,refer.band))
bandsig=unique(bandsig)
genesig=unique(genesig)
LSAres=list()
if (!is.null(bandsig)){
  index=match(as.character(bandsig$cell),as.character(realcell$cell))
  bandsig$subtreesize=realcell$subtreesize[index]
  LSAres$bandLSA=bandsig
  paraBand=table(as.character(bandsig$region))
  paraBand=paraBand[paraBand>1]
  if (length(paraBand)>0){
    paraBandsig=GenePara(bandsig,permuteres,type="band",realcell)
    if (!is.null(paraBandsig)){
      LSAres$paraBand=paraBandsig
    }
  }
}
if (!is.null(genesig)){
  index=match(as.character(genesig$cell),as.character(realcell$cell))
  genesig$subtreesize=realcell$subtreesize[index]
  LSAres$geneLSA=genesig
  paraGene=table(as.character(genesig$region))
  paraGene=paraGene[paraGene>1]
  if (length(paraGene)>0){
      paraGenesig=GenePara(genesig,permuteres,type="gene",realcell)
      LSAres$paraGene=paraGenesig
  }
}
print.noquote("Estimate parallel evolution")
allsig=c()
if ("geneLSA" %in% names(LSAres)){
  geneLSA=LSAres$geneLSA
  geneLSA$CNA="AMP"
  geneLSA$CNA[geneLSA$Score<0]="DEL"
  allsig=rbind(allsig,geneLSA)
  write.table(geneLSA,paste(outpath,"/gene.LSA.txt",sep=""),col.names = T,row.names=F,sep="\t",quote=F)
}else{
  print.noquote("No LSA is identified at gene level!")
}
if ("bandLSA" %in% names(LSAres)){
  bandLSA=LSAres$bandLSA
  bandLSA$CNA="AMP"
  bandLSA$CNA[bandLSA$Score<0]="DEL"
  allsig=rbind(allsig,bandLSA)
  write.table(bandLSA,paste(outpath,"/segmental.LSA.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}else{
  print.noquote("No segmental LSA is identified!")
}
paraEvent=c()
if ("paraBand" %in% names(LSAres)){
  paraEvent=rbind(paraEvent,LSAres$paraBand)
}
if ("paraGene" %in% names(LSAres)){
  paraEvent=rbind(paraEvent,LSAres$paraGene)
}
if ("paraBand" %in% names(LSAres)|"paraGene" %in% names(LSAres)){
  write.table(paraEvent,paste(outpath,"/parallel.LSA.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=FALSE)
}
####plot LSA Tree
if (!is.null(allsig)){
  LSAnetwork=CNAconnect(allsig,celltree)
  nodes=data.frame(id=union(LSAnetwork[,1],LSAnetwork[,2]),size=5)
  tab=table(as.character(allsig$cell))
  index=match(nodes$id,names(tab))
  nodes$size[!is.na(index)]=nodes$size[!is.na(index)]*tab[index[!is.na(index)]]/5
  nodes$size[nodes$size<=5]=5
  nodes$color="gray"
  nodes$color[!is.na(index)]=rainbow(length(unique(allsig$cell)))
  annotation=c()
  for (i in 1:dim(nodes)[1]){
    if (as.character(nodes$id[i]) %in% as.character(allsig$cell)){
      CNA=allsig[as.character(allsig$cell)==as.character(nodes$id[i]),]
      #pvalue=apply(CNA[,4:5],1,min)
      CNA=CNA[order(CNA$pvalue),]
      CNA=paste(as.character(CNA$region),as.character(CNA$CNA),sep=":")
      CNA1=CNA[1]
      if (length(CNA)>1){
        for (j in 2:min(5,length(CNA))){
          CNA1=paste(CNA1,CNA[j],sep=";")
          }
      }
      annotation[i]=CNA1
    }
  }
  nodes$annotation=annotation
  nodes$size=nodes$size/max(nodes$size)*100
  links=data.frame(from=LSAnetwork[,1],to=LSAnetwork[,2],weight=as.numeric(LSAnetwork[,3]))
  pdf(file=paste(outpath,"/LSA.tree.pdf",sep=""),width = 6,height = 6,useDingbats = F)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
  plot(net, layout=layout_as_tree,vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label.cex=0.5,vertex.label=nodes$annotation)
  dev.off()
}
