####permute CNV profile
BinCNV<-function(refgene,newdata,delt){
  index=match(colnames(newdata),as.character(refgene[,1]))
  newgene=data.frame(gene=colnames(newdata),chr=refgene[index,2],start=refgene[index,3],end=refgene[index,4])
  newgene[,2]=as.character(newgene[,2])
  newgene[newgene[,2]=="chrX",2]="chr23"
  chrom=paste("chr",c(1:23),sep="")
  chro=as.character(newgene$chr)
  chrom=intersect(chrom,chro)
  #chrom=c(chrom,"chrX")
  bin=delt
  regionCNV=c()
  chrregion=c()
  for (i in 1:length(chrom)){
    subgene=as.character(newgene$gene[newgene$chr==chrom[i]])
    kk=length(subgene)/bin
    inteKK=round(kk)
    if (inteKK>1){
      for(j in 1:(inteKK-1)){
        subgene1=subgene[((j-1)*bin+1):(j*bin)]
        index=match(subgene1,colnames(newdata))
        regionCNV=rbind(regionCNV,apply(newdata[,index],1,mean))
        chrregion=c(chrregion,paste(chrom[i],"_",j,sep=""))
      }
      subgene1=subgene[((inteKK-1)*bin+1):length(subgene)]
      index=match(subgene1,colnames(newdata))
      regionCNV=rbind(regionCNV,apply(newdata[,index],1,mean))
      chrregion=c(chrregion,paste(chrom[i],"_",inteKK,sep=""))
    }else{
      index=match(subgene,colnames(newdata))
      regionCNV=rbind(regionCNV,apply(newdata[,index],1,mean))
      chrregion=c(chrregion,paste(chrom[i],"_",1,sep=""))
    }
  }
  row.names(regionCNV)=chrregion
  regionCNV=round(regionCNV)
  return(regionCNV)
}

permuteSeg <- function(CNV,generegion){
  CNV=CNV[,3:dim(CNV)[2]]
  row.names(CNV)=generegion$ID1
  CNV=t(CNV)
  chr=colnames(CNV)
  chr=do.call(rbind,strsplit(chr,split="_"))
  chr=unique(chr[,1])
  CNV1=do.call(cbind,lapply(chr,function(j,CNV){
  		subCNV=CNV[,grep(paste(j,"_",sep=""),colnames(CNV))]
  		index=sample(1:dim(subCNV)[1],dim(subCNV)[1])
  		return(subCNV[index,])
  	},CNV=CNV))
  colnames(CNV1)=colnames(CNV)
  row.names(CNV1)=row.names(CNV)
  return(CNV1)
}

permuteGene <- function(CNV,reference){
  CNV=t(CNV)
  index=match(colnames(CNV),as.character(reference[,1]))
  chr=as.character(reference[index,2])
  chr[chr=="chrX"]="chr23"
  chr[chr=="chrY"]="chr24"
  colnames(CNV)=paste(chr,"_",colnames(CNV),sep="")
  chr=colnames(CNV)
  chr=do.call(rbind,strsplit(chr,split="_"))
  chr=unique(chr[,1])
  CNV1=do.call(cbind,lapply(chr,function(j,CNV){
      subCNV=CNV[,grep(paste(j,"_",sep=""),colnames(CNV))]
      index=sample(1:dim(subCNV)[1],dim(subCNV)[1])
      return(subCNV[index,])
  },CNV=CNV))
  genename=do.call(rbind,strsplit(colnames(CNV),split="_"))[,2]
  colnames(CNV1)=genename
  row.names(CNV1)=row.names(CNV)
  return(CNV1)
}

permuteID <- function(ID,permuteCNV,ans){
  permuteCNV1=do.call(cbind,lapply(ID, function(id,cnv,ans){
    chrID=ans$ID1[ans$ID==id]
    index=match(chrID,colnames(cnv))
    if (length(index)==1){
      return(as.numeric(cnv[,index]))
    }else{
      return(round(apply(cnv[,index],1,mean)))
    }
  },cnv=permuteCNV,ans=ans))
  colnames(permuteCNV1)=ID
  return(permuteCNV1)
}

permuteScore <- function(data,ID,ans,datatype,pathwaygene,generegion,reference){
  if (datatype == "D"){
    permuteCNV=permuteSeg(data,generegion)
    permuteCNV1=permuteID(ID,permuteCNV,ans)
    geneCNV=do.call(cbind,lapply(1:dim(pathwaygene)[1],geneCNAfunction,pathwaygene=pathwaygene,ancestorCNV=permuteCNV,generegion=generegion))
    colnames(geneCNV)=as.character(pathwaygene$name)
    rownames(geneCNV)=rownames(permuteCNV)
    index=apply(geneCNV,2,function(x){
      if (NA %in% x){
        return(0)
      }else{
        return(1)
      }
    })
    geneCNV=geneCNV[,index==1]
  }
  if (datatype == "R"){
    permuteCNV=permuteGene(data,reference)
    index=match(colnames(permuteCNV),row.names(generegion))
    permuteCNV1=permuteCNV
    colnames(permuteCNV1)=generegion$ID1[index]
    permuteCNV1=permuteID(ID,permuteCNV1,ans)
    index=match(as.character(pathwaygene$name),colnames(permuteCNV))
    geneCNV=permuteCNV[,index[!is.na(index)]]
  }
  Gscore=lapply(as.character(cell1$cell),lineageScore,permuteCNV1,celltree)
  names(Gscore)=as.character(cell1$cell)
  geneGscore=lapply(as.character(cell1$cell),lineageScore,geneCNV,celltree)
  names(geneGscore)=as.character(cell1$cell)
  return(list(bandGscore=Gscore,geneGscore=geneGscore))
}


#####permutation tree
permuteTreeScore <- function(permutetree,ID,ans,datatype,pathwaygene,generegion,reference,permutationPath){
  j=strsplit(permutetree,split="[.]")[[1]][2]
  if (datatype == "D"){
    permuteCNV=read.csv(paste(permutationPath,"/permute.",j,".CNV.txt",sep=""),sep="\t")
    permuteCNV1=permuteID(ID,permuteCNV,ans)
    geneCNV=do.call(cbind,lapply(1:dim(pathwaygene)[1],geneCNAfunction,pathwaygene=pathwaygene,ancestorCNV=permuteCNV,generegion=generegion))
    colnames(geneCNV)=as.character(pathwaygene$name)
    index=apply(geneCNV,2,function(x){
      if (NA %in% x){
        return(0)
      }else{
        return(1)
      }
    })
    geneCNV=geneCNV[,index==1]
  }
  if (datatype == "R"){
    permuteCNV=read.csv(paste(permutationPath,"/permute.",j,".gene.CNV.txt",sep=""),sep="\t")
    index=match(colnames(permuteCNV),row.names(generegion))
    permuteCNV1=permuteCNV
    colnames(permuteCNV1)=generegion$ID1[index]
    permuteCNV1=permuteID(ID,permuteCNV1,ans)
    index=match(as.character(pathwaygene$name),colnames(permuteCNV))
    geneCNV=permuteCNV[,index[!is.na(index)]]
  }
  permutecelltree=read.csv(paste(permutationPath,"/permute.",j,".CNV.txt.celltree.txt",sep=""),sep="\t")
  permutecell=union(as.character(permutecelltree[,1]),as.character(permutecelltree[,2]))
  permutecell=data.frame(cell=permutecell)
  permutecell$depth=sapply(as.character(permutecell$cell),depthFunction,cellTree=permutecelltree)
  permutecell$height=sapply(as.character(permutecell$cell),heightFunction,cellTree=permutecelltree)
  permutecell$subtreesize=sapply(as.character(permutecell$cell),subtreeSize,cellTree=permutecelltree)
  permutecell1=permutecell[permutecell$subtreesize>=5,]
  permutecell1=permutecell1[permutecell1$cell!="root",]
  Gscore=lapply(as.character(permutecell1$cell),lineageScore,permuteCNV1,permutecelltree)
  names(Gscore)=as.character(permutecell1$cell)
  geneGscore=lapply(as.character(permutecell1$cell),lineageScore,geneCNV,permutecelltree)
  names(geneGscore)=as.character(permutecell1$cell)
  return(list(cell=permutecell1,bandGscore=Gscore,geneGscore=geneGscore))
}

PermuteSig <- function(observe,background,type){
  index=match(names(observe),colnames(background))
  observe=observe[!is.na(index)]
  background=background[,index[!is.na(index)]]
  if (type=="Amp"){
    pvalue=unlist(lapply(1:length(observe), function(j,observe,background){
      if (is.null(dim(background))){
        x=background[j]
        p=(length(x[x>=observe[j]])+1)/(length(x)+1)
      }else{
        x=background[,j]
        p=(length(x[x>=observe[j]])+1)/(length(x)+1)
      }
      return(p)
    },observe,background))
  }
  if (type=="Del"){
    pvalue=unlist(lapply(1:length(observe), function(j,observe,background){
      if (is.null(dim(background))){
        x=background[j]
        p=(length(x[x>=observe[j]])+1)/(length(x)+1)
      }else{
        x=background[,j]
        p=(length(x[x>=observe[j]])+1)/(length(x)+1)
      }
      return(p)
    },abs(observe),abs(background)))
  }
  names(pvalue)=names(observe)
  return(pvalue)
}


significanceLevel <- function(x,realband,realgene,permuteres,realcell){
  sizecutoff=realcell$subtreesize[realcell$cell==as.character(realcell$cell)[x]]
  if ("cell" %in% names(permuteres[[1]])){
    permutedis=lapply(1:length(permuteres),function(j,permuteres,cutoff){
      permutecell=permuteres[[j]]$cell
      delt=abs(permutecell$subtreesize-cutoff)
      permuteindex=which(delt==min(delt))
      permutebandA=c()
      permutebandD=c()
      permutegeneA=c()
      permutegeneD=c()
      for (i in permuteindex){
        permutebandA=rbind(permutebandA,permuteres[[j]]$bandGscore[[i]][1,])
        permutebandD=rbind(permutebandD,permuteres[[j]]$bandGscore[[i]][2,])
        permutegeneA=rbind(permutegeneA,permuteres[[j]]$geneGscore[[i]][1,])
        permutegeneD=rbind(permutegeneD,permuteres[[j]]$geneGscore[[i]][2,])
      }
      return(list(permutebandA,permutebandD,permutegeneA,permutegeneD))
    },permuteres,sizecutoff)
  }else{
    permutedis=lapply(1:length(permuteres),function(j,permuteres,x){
        permutebandA=permuteres[[j]]$bandGscore[[x]][1,]
        permutebandD=permuteres[[j]]$bandGscore[[x]][2,]
        permutegeneA=permuteres[[j]]$geneGscore[[x]][1,]
        permutegeneD=permuteres[[j]]$geneGscore[[x]][2,]
        return(list(permutebandA,permutebandD,permutegeneA,permutegeneD))
    },permuteres,x)
  }
  bandAdis=do.call(rbind,lapply(1:length(permutedis), function(j,permutedis){
    return(permutedis[[j]][[1]])
  },permutedis))
  bandDdis=do.call(rbind,lapply(1:length(permutedis), function(j,permutedis){
    return(permutedis[[j]][[2]])
  },permutedis))
  geneAdis=do.call(rbind,lapply(1:length(permutedis), function(j,permutedis){
    return(permutedis[[j]][[3]])
  },permutedis))
  geneDdis=do.call(rbind,lapply(1:length(permutedis), function(j,permutedis){
    return(permutedis[[j]][[4]])
  },permutedis))
  pvalueAband=PermuteSig(realband[[x]][1,],bandAdis,type="Amp")
  pvalueDband=PermuteSig(realband[[x]][2,],bandDdis,type="Del")
  pvalueAgene=PermuteSig(realgene[[x]][1,],geneAdis,type="Amp")
  pvalueDgene=PermuteSig(realgene[[x]][2,],geneDdis,type="Del")
  adjustAband=p.adjust(pvalueAband, method = "BH")
  adjustDband=p.adjust(pvalueDband, method = "BH")
  adjustAgene=p.adjust(pvalueAgene, method = "BH")
  adjustDgene=p.adjust(pvalueDgene, method = "BH")
  index=match(names(realband[[x]][1,]),names(pvalueAband))
  bandAscore=realband[[x]][1,][!is.na(index)]
  bandDscore=realband[[x]][2,][!is.na(index)]
  index=match(names(realgene[[x]][1,]),names(pvalueAgene))
  geneAscore=realgene[[x]][1,][!is.na(index)]
  geneDscore=realgene[[x]][2,][!is.na(index)]
  Aband=data.frame(region=names(bandAscore),Score=bandAscore,pvalue=pvalueAband,adjustp=adjustAband)
  Dband=data.frame(region=names(bandDscore),Score=bandDscore,pvalue=pvalueDband,adjustp=adjustDband)
  Agene=data.frame(region=names(geneAscore),Score=geneAscore,pvalue=pvalueAgene,adjustp=adjustAgene)
  Dgene=data.frame(region=names(geneDscore),Score=geneDscore,pvalue=pvalueDgene,adjustp=adjustDgene)
  return(list(Aband,Dband,Agene,Dgene))
}


RefineCNA <- function(res,ancestor,celltree,realcell){
  res1=c()
  res$ID=paste(res$cell,res$region,res$Score,sep=":")
  for (i in 1:length(ancestor)){
    region=as.character(res$region[as.character(res$cell)==ancestor[i]])
    if (length(region)!=0){
      child=splitTree(ancestor[i],celltree = celltree)
      overlap=intersect(setdiff(child,ancestor[i]),as.character(res$cell))
      if (length(overlap)==0){
        res1=rbind(res1,res[res$cell==ancestor[i],])
      }else{
        subcell=c(ancestor[i],overlap)
        filterregion=lapply(region, function(r,res,subcell,realcell,child){
          subres=res[res$region==r,]
          index=match(as.character(subres$cell),subcell)
          subres=subres[!is.na(index),]
          if (dim(subres)[1]==1){
            anchorIndex=1
          }else{
            index=match(child,realcell$cell)
            subchild=realcell[index[!is.na(index)],]
            index=match(subchild$cell,subres$cell)
            if (sum(is.na(index))==0){
              anchorIndex=which(subres$cell==ancestor[i])
            }else{
              M=min(subchild$depth[is.na(index)])
              if(length(which(subres$depth>=M))>0){
                anchorIndex=min(which(subres$depth>=M))
              }else{
                anchorIndex=which(subres$cell==ancestor[i])
              }
            }
          }
          return(list(keep=subres[anchorIndex,],remove=subres[-anchorIndex,]))
        },res,subcell,realcell,child)
        keepR=do.call(rbind,lapply(1:length(filterregion), function(j,filterregion){return(filterregion[[j]]$keep)},filterregion))
        removeR=do.call(rbind,lapply(1:length(filterregion), function(j,filterregion){return(filterregion[[j]]$remove)},filterregion))
        keepR=keepR[!is.na(keepR$region),]
        removeR=removeR[!is.na(removeR$region),]
        index=match(removeR$ID,res$ID)
        if (length(index)!=0){
          res=res[-index,]
        }
        index=match(keepR$ID,res$ID)
        if (length(index)!=0){
          res=res[-index,]
        }
        res1=rbind(res1,keepR)
      }
    }
  }
  res1=res1[order(res1$depth),1:6]
  res1=MergeLSA(res1,celltree,realcell)
  res1=res1[order(res1$depth),1:6]
  return(res1)
}

MergeLSA <- function(res1,celltree,realcell){
  region=table(as.character(res1$region))
  region1=names(region[region==1])
  region2=names(region[region>1])
  fres=c()
  if (length(region1)>0){
    index=match(region1,as.character(res1$region))
    fres=rbind(fres,res1[index,])
  }
  if (length(region2)>0){
    for (i in 1:length(region2)){
      subres=res1[res1$region==region2[i],]
      while(1){
        child=unlist(lapply(as.character(subres$cell), function(x,subres,celltree){
          return(setdiff(splitTree(x,celltree),x))
        },subres,celltree))
        index=match(intersect(child,as.character(subres$cell)),as.character(subres$cell))
        if (length(index)>0){
          subres=subres[-index,]
        }
        index=match(as.character(subres$cell),celltree[,2])
        parent=as.character(celltree[index,1])
        parentNum=table(parent)
        parentNum=parentNum[parentNum>1]
        if ("root" %in% names(parentNum)){
          parentNum=parentNum[names(parentNum)!="root"]
        }
        if (length(parentNum)>0){
          ancestor=names(parentNum)
          addres=do.call(rbind,lapply(ancestor, function(a,subres,parent,region){
            kk=which(parent==a)
            addres=apply(subres[kk,2:4],2,mean)
            addres=c(region2[i],addres,a,mean(subres$depth[kk])-1)
            return(addres)
          },subres,parent,region))
          colnames(addres)=colnames(subres)
          index=match(parent,ancestor)
          kk=which(!is.na(index))
          subres=subres[-kk,]
          subres=rbind(subres,addres)
          subres$Score=as.numeric(as.character(subres$Score))
          subres$pvalue=as.numeric(as.character(subres$pvalue))
          subres$adjustp=as.numeric(as.character(subres$adjustp))
          subres$depth=as.numeric(as.character(subres$depth))
        }else{
          break
        }
      }
      index=match(as.character(subres$cell),celltree[,2])
      parent=as.character(celltree[index,1])
      index=match(parent,celltree[,2])
      updateparent=parent
      updateparent[is.na(index)]="root"
      if (length(unique(updateparent))==1){
        if (unique(updateparent)=="root"){
          child=as.character(celltree[celltree[,1]==parent[is.na(index)],2])
          overlap=intersect(child,realcell$cell)
          if (length(overlap)>dim(subres)[1]){
            fres=rbind(fres,subres)
          }
        }else{
          fres=rbind(fres,subres)
        }
      }else{
        fres=rbind(fres,subres)
      }
    }
  }
  return(fres)
}

CollectAsso <- function(significance,cutoff,celltree,realcell){
  bandsig=c()
  genesig=c()
  for (i in 1:length(significance)){
    bandres=c()
    generes=c()
    bandres=rbind(significance[[i]][[1]][significance[[i]][[1]][,3]<cutoff,],significance[[i]][[2]][significance[[i]][[2]][,3]<cutoff,])
    bandres=bandres[abs(bandres$Score)>0.5,]
    generes=rbind(significance[[i]][[3]][significance[[i]][[3]][,3]<cutoff,],significance[[i]][[4]][significance[[i]][[4]][,3]<cutoff,])
    generes=generes[abs(generes$Score)>0.5,]
    if (dim(bandres)[1]>0){
      bandres$cell=as.character(realcell$cell[i])
      bandres$depth=realcell$depth[i]
    }
    if (dim(generes)[1]>0){
      generes$cell=as.character(realcell$cell[i])
      generes$depth=realcell$depth[i]
    }
    bandsig=rbind(bandsig,bandres)
    genesig=rbind(genesig,generes)
  }
  bandsig=bandsig[!is.na(bandsig$region),]
  genesig=genesig[!is.na(genesig$region),]
  bandsig1=c()
  genesig1=c()
  if (dim(bandsig)[1]>0){
    bandsig=unique(bandsig[order(bandsig$depth),])
    ancestor=unique(as.character(bandsig$cell))
    if (dim(bandsig[bandsig$Score>0,])[1]>0){
      bandsig1=rbind(bandsig1,RefineCNA(bandsig[bandsig$Score>0,],unique(as.character(bandsig$cell[bandsig$Score>0])),celltree,realcell))
    }
    if (dim(bandsig[bandsig$Score<0,])[1]>0){
      bandsig1=rbind(bandsig1,RefineCNA(bandsig[bandsig$Score<0,],unique(as.character(bandsig$cell[bandsig$Score<0])),celltree,realcell))
    }
    bandsig1=bandsig1[order(bandsig1$depth),]
  }
  if (dim(genesig)[1]>0){
    genesig=unique(genesig[order(genesig$depth),])
    ancestor=unique(as.character(genesig$cell))
    if (dim(genesig[genesig$Score>0,])[1]>0){
      genesig1=rbind(genesig1,RefineCNA(genesig[genesig$Score>0,],unique(as.character(genesig$cell[genesig$Score>0])),celltree,realcell))
    }
    if (dim(genesig[genesig$Score<0,])[1]>0){
      genesig1=rbind(genesig1,RefineCNA(genesig[genesig$Score<0,],unique(as.character(genesig$cell[genesig$Score<0])),celltree,realcell))
    }
    genesig1=genesig1[order(genesig1$depth),]
  }
  bandsig1=bandsig1[!is.na(bandsig1$Score),]
  genesig1=genesig1[!is.na(genesig1$Score),]
  return(list(bandres=unique(bandsig1),generes=unique(genesig1)))
}

GenePara <- function(sigres,permuteres,type,realcell){
  paraGene=table(as.character(sigres$region))
  paraGene=paraGene[paraGene>1]
  pvalue=unlist(lapply(1:length(paraGene),function(x,paraGene,sigres){
    region=names(paraGene)[x]
    subres=sigres[sigres$region==region,]
    CNAsign=sign(subres$Score)
    if (length(CNAsign[CNAsign==1])>1){
      subsize=subres$subtreesize[CNAsign==1]
      CNAindex=1
    }
    if (length(CNAsign[CNAsign==(-1)])>1){
      subsize=subres$subtreesize[CNAsign==(-1)]
      CNAindex=2
    }
    if (length(CNAsign[CNAsign==1])>1 | length(CNAsign[CNAsign==(-1)])>1){
      permutelineage=lapply(1:length(permuteres), function(i,permuteres,subsize,realcell,CNAindex){
        if ("cell" %in% names(permuteres[[i]])){
          permutecell=permuteres[[i]]$cell
        }else{
          permutecell=realcell
        }
        if (type=="gene"){
          permuteScore=permuteres[[i]]$geneGscore
        }
        if (type=="band"){
          permuteScore=permuteres[[i]]$bandGscore
        }
        permutelineage=c()
        for (s in 1:length(subsize)){
          delt=abs(permutecell$subtreesize-subsize[s])
          MM=which(delt==min(delt))
          score1=c()
          if (CNAindex == 1){
            for (j in MM){
              score1=c(score1,max(permuteScore[[j]][CNAindex,grep(region,colnames(permuteScore[[j]]))]))
            }
            permutelineage = c(permutelineage,length(score1[score1>=0.5]))
          }
          if (CNAindex == 2){
            for (j in MM){
              score1=c(score1,min(permuteScore[[j]][CNAindex,grep(region,colnames(permuteScore[[j]]))]))
            }
            permutelineage=c(permutelineage,length(score1[score1<=(-0.5)]))
          }
        }
        return(permutelineage)
      },permuteres,subsize,realcell,CNAindex)
      permutelineage=sapply(permutelineage,function(s){return(sum(s))})
      pvalue=(length(permutelineage[permutelineage>=length(subsize)])+1)/(length(permutelineage)+1)
      #pvalue=ppois(length(subsize), mean(permutelineage), lower.tail = F, log.p = FALSE)
      return(pvalue)
    }else{
      return(NA)
    }
  },paraGene,sigres))
  paraGene=data.frame(region=names(paraGene),lineage=as.numeric(paraGene),pvalue=pvalue)
  return(paraGene)
}

mergeCNA <- function(node,sigCNA,band.region,arm.band,refer.band){
  refineSigregion <- function(i,psub,type,refer.band,arm.band,subres,subregion){
    if (dim(psub)[1]==1){
      subres=rbind(subres,psub)
    }else if (dim(psub)[1]>1){
      index=match(as.character(psub$region),as.character(refer.band$ID))
      subrefer=refer.band[index,]
      sublength=sum(subrefer$V3-subrefer$V2)
      f=sublength/arm.band$length[arm.band$chr==i&arm.band$band=="p"]
      if (f > 0.5){
        subres=rbind(subres,data.frame(region=paste(i,type,sep=":"),Score=mean(psub$Score),pvalue=mean(psub$pvalue),adjustp=mean(psub$adjustp),cell=unique(as.character(psub$cell)),depth=unique(psub$depth)))
      }else{
        subband=subregion[subregion[,1]==i,]
        subband=do.call(rbind,strsplit(subband[,2],split="[.]"))
        subband=unique(subband[,1])
        sub1=do.call(rbind,lapply(subband,function(x,psub){
          sub1=psub[grep(x,as.character(psub$region)),]
          if (dim(sub1)[1]==1){
            return(sub1)
          }else{
            index=match(as.character(sub1$region),as.character(refer.band$ID))
            f=sum(refer.band$V3[index]-refer.band$V2[index])/band.region$length[band.region$chr==i&band.region$band==x]
            if (f > 0.5){
              return(data.frame(region=paste(i,x,sep=":"),Score=mean(sub1$Score),pvalue=mean(sub1$pvalue),adjustp=mean(sub1$adjustp),cell=unique(as.character(sub1$cell)),depth=unique(sub1$depth)))
            }else{
              return(sub1)
            }
          }
        },psub))
        subres=rbind(subres,sub1)
      }
    }
    return(subres)
  }
  subsig=sigCNA[as.character(sigCNA$cell)==node,]
  subregion=do.call(rbind,strsplit(as.character(subsig$region),split=":"))
  chr=unique(as.character(subregion[,1]))
  subres=c()
  for (i in chr){
    subdata=subsig[subregion[,1]==i,]
    psub=subdata[grep("p",subdata$region),]
    qsub=subdata[grep("q",subdata$region),]
    if (dim(psub)[1]>0){
      DEL=psub[psub$Score<0,]
      AMP=psub[psub$Score>0,]
      if (dim(DEL)[1]>0){
        subres=refineSigregion(i,DEL,type="p",refer.band,arm.band,subres,subregion)
      }
      if (dim(AMP)[1]>0){
        subres=refineSigregion(i,AMP,type="p",refer.band,arm.band,subres,subregion)
      }
    }
    if (dim(qsub)[1]>0){
      DEL=qsub[qsub$Score<0,]
      AMP=qsub[qsub$Score>0,]
      if (dim(DEL)[1]>0){
        subres=refineSigregion(i,DEL,type="q",refer.band,arm.band,subres,subregion)
      }
      if (dim(AMP)[1]>0){
        subres=refineSigregion(i,AMP,type="q",refer.band,arm.band,subres,subregion)
      }
    }
  }
  return(subres)
}

depthFunction <- function(cell,cellTree){
  root=setdiff(as.character(cellTree[,1]),as.character(cellTree[,2]))
  if (cell == root){
    return(0)
  }else{
    k=1
    while(1){
      cell1=as.character(cellTree[as.character(cellTree[,2])==cell,1])
      if (cell1==root){
        break
      }else{
        k=k+1
        cell=cell1
      }
    }
    return(k)
  }
}

heightFunction <- function(cell,cellTree){
  leaf=setdiff(as.character(cellTree[,2]),as.character(cellTree[,1]))
  if (cell %in% leaf){
    return(0)
  }else{
    k=1
    while (1){
      cell1=as.character(cellTree[!is.na(match(as.character(cellTree[,1]),cell)),2])
      index=match(cell1,leaf)
      if (length(index[is.na(index)])==0){
        break
      }else{
        cell=cell1[is.na(index)]
        k=k+1
      }
    }
    return(k)
  }
}


###calculated subtree size for each node
subtreeSize <- function(cell,cellTree){
  leaf=setdiff(as.character(cellTree[,2]),as.character(cellTree[,1]))
  if (cell %in% leaf){
    return(1)
  }else{
    k=1
    while(1){
      cell1=as.character(cellTree[!is.na(match(as.character(cellTree[,1]),cell)),2])
      k=k+length(cell1)
      cell1=setdiff(cell1,leaf)
      if (length(cell1)==0){
        break
      }else{
        cell=cell1
      }
    }
    return(k)
  }
}

#####Searching subtree rooted at node
splitTree <- function(node,celltree){
  child=node
  while(1){
    index=match(as.character(celltree[,1]),node)
    if (sum(!is.na(index))==0){
      break
    }else{
      node=as.character(celltree[!is.na(index),2])
      child=union(child,node)
    }
  }
  return(child)
}

#######Calculate CFL score
lineageScore <- function(node,newCNV,celltree){
  child=splitTree(node,celltree)
  index=match(child,row.names(newCNV))
  subcnv=newCNV[index,]
  Gscore=apply(subcnv,2,function(x){
    f1=length(x[x>2])
    f2=length(x[x<2])
    if (f1>1){
      Gamp=f1*mean(x[x>2]-2)/length(x)
    }else if (f1==1){
      Gamp=f1*(x[x>2]-2)/length(x)
    }else{
      Gamp=f1
    }
    if (f2>1){
      Gdel=f2*mean(x[x<2]-2)/length(x)
    }else if (f2 == 1){
      Gdel=f2*(x[x<2]-2)/length(x)
    }else{
      Gdel=f2
    }
    return(c(Gamp,Gdel))
  })
  return(Gscore)
}

###Gene level copy number profile
geneCNAfunction<-function(k,pathwaygene,ancestorCNV,generegion){
  gene=pathwaygene[k,]
  chr=as.character(gene[1,1])
  start=gene[1,2]
  end=gene[1,3]
  gene=as.character(gene[1,4])
  subregion=generegion[as.character(generegion$chrom)==chr,]
  if (dim(subregion)[1]>0){
    subCNV=ancestorCNV[,grep(paste(chr,"_",sep=""),colnames(ancestorCNV))]
    startDis=abs(subregion$chrompos-start)
    endDis=abs(subregion$chrompos-end)
    startIndex=which.min(startDis)
    endIndex=which.min(endDis)
    if (startIndex==endIndex){
      return(subCNV[,startIndex])
    }else{
      m=min(startIndex,endIndex)
      M=max(startIndex,endIndex)
      subCNV1=subCNV[,m:M]
      a=apply(subCNV1,1,function(x){
        if (length(x[x!=2])==0){
          return(2)
        }else{
          if (length(x[x>2])>length(x[x<2])){
            return(max(x[x>2]))
          }else if (length(x[x>2])<length(x[x<2])){
            return(min(x[x<2]))
          }else{
            if (startDis[startIndex]<endDis[endIndex]){
              return(x[m])
            }else{
              return(x[M])
            }
          }
        }
      })
      return(a)
    }
  }else{
    return(rep(NA,length=dim(ancestorCNV)[1]))
  }
}
####

CNAconnect <- function(sigCNA,celltree){
  node=unique(sigCNA[,5:6])
  node=node[order(node$depth,decreasing = TRUE),]
  CNAnetwork=c()
  for (i in 1:dim(node)[1]){
    if (i == dim(node)[1]){
      node1="root"
      dist=as.numeric(as.character(node$depth[i]))
      CNAnetwork=rbind(CNAnetwork,c(node1,as.character(node$cell[i]),dist))
    }else{
      ancestor=c()
      node2=as.character(node$cell[i])
      while(1){
        node1=as.character(celltree$from[celltree$to==node2])
        if (length(node1)==0){
          break
        }else{
          ancestor=union(ancestor,node1)
          node2=node1
        }
      }
      index=match(as.character(node$cell),ancestor)
      k=which.max(as.numeric(as.character(node$depth[!is.na(index)])))
      if (length(k)==0){
        node1="root"
        dist=as.numeric(as.character(node$depth[i]))
      }else{
        node1=as.character(node[!is.na(index),][1,1])
        dist=as.numeric(as.character(node[i,2]))-as.numeric(as.character(node[!is.na(index),][1,2]))
      }
      CNAnetwork=rbind(CNAnetwork,c(node1,as.character(node$cell[i]),dist))
    }
  }
  colnames(CNAnetwork)=c("node1","node2","dist")
  return(CNAnetwork)
}

CombineRegion <- function(node,newsig,refer.band){
  cellsig=c()
  subsig=newsig[newsig$cell==node,]
  indexregion=match(subsig$region,refer.band$ID)
  if (length(indexregion[!is.na(indexregion)])>0){
    subsig1=subsig[!is.na(indexregion),]
    if (dim(subsig1)[1]==1){
      cellsig=rbind(cellsig,subsig1)
    }else if (dim(subsig1)[1]>1){
      subsig1$band=indexregion[!is.na(indexregion)]
      subsig1=subsig1[order(subsig1$band),]
      subregion=do.call(rbind,sapply(as.character(subsig1$region),strsplit,split=":"))
      subsig1$chr=subregion[,1]
      subsig1$arm=sapply(as.character(subregion[,2]),substr,start=1,stop=1)
      subsig1$ll=sapply(as.character(subregion[,2]),nchar)
      subsig1$bandID=substr(as.character(subregion[,2]),start=2,stop=subsig1$ll)
      k=1
      newsig11=c()
      start=k
      while (k < dim(subsig1)[1]){
        if (subsig1$band[k+1]-subsig1$band[k]==1&(subsig1$Score[k+1]*subsig1$Score[k])>0&subsig1$arm[k+1]==subsig1$arm[k]){
          k=k+1
          if (k == dim(subsig1)[1]){
            end=k
            region=paste(subsig1$chr[start],":",subsig1$arm[start],subsig1$bandID[start],"-",subsig1$bandID[end],sep="")
            tempsig=data.frame(region=region,Score=mean(subsig1$Score[start:end]),pvalue=mean(subsig1$pvalue[start:end]),adjustp=mean(subsig1$adjustp[start:end]),cell=subsig1$cell[start],depth=subsig1$depth[start],subtreesize=subsig1$subtreesize[start])
            newsig11=rbind(newsig11,tempsig)
          }
        }else{
          end=k
          if (start == end){
            newsig11=rbind(newsig11,subsig1[k,1:7])
          }else{
            region=paste(subsig1$chr[start],":",subsig1$arm[start],subsig1$bandID[start],"-",subsig1$bandID[end],sep="")
            tempsig=data.frame(region=region,Score=mean(subsig1$Score[start:end]),pvalue=mean(subsig1$pvalue[start:end]),adjustp=mean(subsig1$adjustp[start:end]),cell=subsig1$cell[start],depth=subsig1$depth[start],subtreesize=subsig1$subtreesize[start])
            newsig11=rbind(newsig11,tempsig)
          }
          k=k+1
          start=k
          if (start==dim(subsig1)[1]){
            newsig11=rbind(newsig11,subsig1[k,1:7])
          }
        }
      }
      cellsig=rbind(cellsig,newsig11)
    }
  }
  if (length(indexregion[is.na(indexregion)])>0){
    subsig2=subsig[is.na(indexregion),]
    if (dim(subsig2)[1]==1){
      cellsig=rbind(cellsig,subsig2[,1:7])
    }else if (dim(subsig2)[1]>1){
      subregion2=do.call(rbind,sapply(as.character(subsig2$region),strsplit,split=":"))
      index=match(subregion2[,2],c("p","q"))
      if (length(index[!is.na(index)])>0){
        subsig21=subsig2[!is.na(index),]
        cellsig=rbind(cellsig,subsig21[,1:7])
      }
      if (length(index[is.na(index)])>0){
        subsig22=subsig2[is.na(index),]
        subsig22=subsig22[order(as.character(subsig22$region)),]
        subregion=do.call(rbind,sapply(as.character(subsig22$region),strsplit,split=":"))
        subsig22$chr=subregion[,1]
        subsig22$arm=sapply(as.character(subregion[,2]),substr,start=1,stop=1)
        subsig22$ll=sapply(as.character(subregion[,2]),nchar)
        subsig22$band=as.numeric(sapply(as.character(subregion[,2]),substr,start=2,stop=subsig22$ll))
        k=1
        newsig22=c()
        start=k
        while (k < dim(subsig22)[1]){
          if (subsig22$chr[k+1]==subsig22$chr[k]&subsig22$arm[k+1]==subsig22$arm[k]&(subsig22$Score[k+1]*subsig22$Score[k])>0&subsig22$band[k+1]-subsig22$band[k]==1){
            k=k+1
            if (k == dim(subsig22)[1]){
              end=k
              region=paste(subsig22$chr[start],":",subsig22$arm[start],subsig22$band[start],"-",subsig22$band[end],sep="")
              tempsig=data.frame(region=region,Score=mean(subsig22$Score[start:end]),pvalue=mean(subsig22$pvalue[start:end]),adjustp=subsig22$adjustp[start:end],cell=subsig22$cell[start],depth=subsig22$depth[start],subtreesize=subsig22$subtreesize[start])
              newsig22=rbind(newsig22,tempsig)
            }
          }else{
            end=k
            if (start == end){
              newsig22=rbind(newsig22,subsig22[k,1:7])
            }else{
              region=paste(subsig22$chr[start],":",subsig22$arm[start],subsig22$band[start],"-",subsig22$band[end],sep="")
              tempsig=data.frame(region=region,Score=mean(subsig22$Score[start:end]),pvalue=mean(subsig22$pvalue[start:end]),adjustp=subsig22$adjustp[start:end],cell=subsig22$cell[start],depth=subsig22$depth[start],subtreesize=subsig22$subtreesize[start])
              newsig22=rbind(newsig22,tempsig)
            }
            k=k+1
            start=k
            if (start==dim(subsig22)[1]){
              newsig22=rbind(newsig22,subsig22[k,1:7])
            }
          }
        }
        cellsig=rbind(cellsig,newsig22)
      }
    }
  }
  return(cellsig)
}
