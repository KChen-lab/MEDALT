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

#######Association test
MonotoneTest <- function(j,data){
  model1=wilcox.test(as.numeric(data[data$cohort==0,j]),as.numeric(data[data$cohort==1,j]))
  return(model1$p.value)
}

TrendTest <- function(j,data){
  tab=table(data[,j], data$cohort)
  model2=CochranArmitageTest(tab)
  return(model2$p.value)
}

LikelihoodRatio <- function(j, data){
  signal=data[,j]
  batch=factor(data$cohort)
  sample=factor(data$subject)
  trait <- data$cohort
  set.seed(0)
  results <- CNVtest.select.model(signal=signal, batch = batch, sample = sample, n.H0 = 3, method="BIC", v.ncomp = 1:5, v.model.component = rep('gaussian',5), v.model.mean = rep("~ strata(cn)",5), v.model.var = rep("~1", 5))
  ncomp <- results$selected
  if (ncomp == 1){
    ncomp = 2
  }
  fit.pca <- CNVtest.binary ( signal = signal, sample = sample, batch = batch, ncomp = ncomp, n.H0=3, n.H1=0, model.var= '~ strata(cn)')
  pca.posterior <- as.matrix((fit.pca$posterior.H0)[, paste('P',seq(1:ncomp),sep='')])
  dimnames(pca.posterior)[[1]] <- (fit.pca$posterior.H0)$subject
  signal=data.frame(signal=signal)
  signal=as.matrix(signal)
  dimnames(signal)[[1]]=data$subject
  ldf.signal <- apply.ldf(signal, pca.posterior)
  fit.ldf <- try(CNVtest.binary ( signal = signal, sample = sample, batch = batch, disease.status = trait, ncomp = ncomp, n.H0=3, n.H1=1, model.var = "~cn"),silent=TRUE)
  if ("model.H0" %in% names(fit.ldf)){
    LR.statistic <- -2*(fit.ldf$model.H0$lnL - fit.ldf$model.H1$lnL)
    return(pchisq(LR.statistic,df=1,lower.tail = F))
  }else{
    return(NA)
  }
}
####
CNAassociation <- function(cnv,cell,celltree,method="Standard"){
  index=apply(cnv,2,function(x){
    y=table(x)
    if (length(y)==1){
      return(0)
    }else{
      if (min(y)<5&length(y)==2){
        return(0)
      }else{
        return(1)
      }
    }
  })
  cnv1=cnv[,index==1]
  subroot=cell[cell$subtreesize>=5,]
  if (dim(subroot)[1]>0){
    depth=unique(subroot$depth)
    depth=depth[order(depth)]
    depth=setdiff(depth,0)
    AssociateRes<-lapply(depth,function(d,subroot,celltree,cnv1,method){
      ancestor=as.character(subroot$cell[subroot$depth==d])
      child=lapply(ancestor,splitTree,celltree=celltree)
      pvalue=lapply(1:length(ancestor), function(i,child,cell,cnv1,method){
        case=child[[i]]
        control=setdiff(as.character(cell$cell[cell$depth!=0]),case)
        if (length(control)>=5){
          #ancestorIndex[i]=1
          data=data.frame(subject=c(case,control),cohort=c(rep(1,length=length(case)),rep(0,length=length(control))))
          index=match(c(case,control),row.names(cnv1))
          data=cbind(data,cnv1[index,])
          model1res=sapply(3:dim(data)[2], MonotoneTest,data=data)
          model2res=sapply(3:dim(data)[2], TrendTest,data=data)
          #model1.pvalue=model1res
          #model2.pvalue=model2res
          model1.pvalue=p.adjust(model1res)
          model2.pvalue=p.adjust(model2res)
          if (method=="Standard"){
            pvalue=data.frame(region=colnames(cnv1),monotone.test=model1.pvalue,trend.test=model2.pvalue)
          }else{
            model3res=sapply(3:dim(data)[2], LikelihoodRatio,data=data)
            model3.pvalue=p.adjust(model3res)
            pvalue=data.frame(region=colnames(cnv1),monotone.test=model1.pvalue,trend.test=model2.pvalue,likelihood.test=model3.pvalue)
          }
          return(pvalue)
        }
      },child=child,cell=cell,cnv1=cnv1,method=method)
      ancestorIndex=sapply(1:length(ancestor), function(i,child,cell){
        case=child[[i]]
        control=setdiff(as.character(cell$cell[cell$depth!=0]),case)
        if (length(control)!=0){
          return(1)
        }else{
          return(0)
        }
      },child=child,cell=cell)
      if (length(ancestorIndex)!=1){
        names(pvalue)=ancestor[ancestorIndex==1]
      }else{
        names(pvalue)=ancestor
      }
      #names(pvalue)=ancestor[ancestorIndex==1]
      return(pvalue)
    },subroot=subroot,celltree=celltree,cnv1=cnv1,method=method)
    names(AssociateRes)=depth
    return(AssociateRes)
  }
}

CollectAsso <- function(pvalue,cutoff,cnv,celltree){
  res=c()
  for (i in 1:length(pvalue)){
    for (j in 1:length(pvalue[[i]])){
      if (!is.null(pvalue[[i]][[j]])){
        minP=apply(pvalue[[i]][[j]][,2:3],1,min)
        if (sum(minP<cutoff)>0){
          res=rbind(res,cbind(names(pvalue)[i],names(pvalue[[i]])[j],pvalue[[i]][[j]][minP<cutoff,]))
        }
      }
    }
  }
  if (!is.null(res)){
    colnames(res)=c("depth","cell",colnames(res)[3:dim(res)[2]])
    ancestor=unique(as.character(res$cell))
    res1=c()
    for (i in 1:length(ancestor)){
      region=as.character(res$region[as.character(res$cell)==ancestor[i]])
      if (length(region)!=0){
        child=splitTree(ancestor[i],celltree = celltree)
        overlap=intersect(setdiff(child,ancestor[i]),as.character(res$cell))
        if (length(overlap)==0){
          res1=rbind(res1,res[res$cell==ancestor[i],])
        }else{
          subres=do.call(rbind,lapply(c(ancestor[i],overlap),function(node,cnv,celltree){
            subchi=splitTree(node,celltree = celltree)
            index1=match(subchi,row.names(cnv))
            index2=match(region,colnames(cnv))
            subcnv=cnv[index1,index2]
            if (length(index2)==1){
              fre=length(subcnv[subcnv!=2])/length(subcnv)
            }else{
              fre=apply(subcnv,2,function(x){
                return(length(x[x!=2])/length(x))
              })
            }
            return(fre)
          },cnv=cnv,celltree=celltree))
          row.names(subres)=c(ancestor[i],overlap)
          ancestral=apply(subres,2,function(x){
            MM=max(x)
            index=which(x==MM)
            return(c(ancestor[i],overlap)[index])
          })
          for (j in 1:length(ancestral)){
            if (length(ancestral[[j]])==1){
              node1=ancestral[[j]]
              res1=rbind(res1,res[res$cell==node1&res$region==names(ancestral)[j],])
            }else{
              index=match(ancestral[[j]],as.character(cell$cell))
              node1=ancestral[[j]][which.min(cell$depth[index])]
              res1=rbind(res1,res[res$cell==node1&res$region==names(ancestral)[j],])
            }
            index1=match(as.character(res$cell),setdiff(c(ancestor[i],overlap),node1))
            index2=match(as.character(res$region),names(ancestral)[j])
            if (length(intersect(which(!is.na(index1)),which(!is.na(index2))))!=0){
              res=res[-intersect(which(!is.na(index1)),which(!is.na(index2))),]
            }
          }
        }

      }

    }
    res1=res1[order(res1$depth),]
    res1=unique(res1)
    return(res1)
  }
}

RefineCNA <- function(cnv,cell,celltree,method="Standard",cutoff=0.01){
  pvalue=CNAassociation(cnv,cell,celltree)
  #pvalue=CNAassociation(cnv=cnv,cell=cell,celltree = celltree,method = "all")
  if (length(pvalue)!=0){
    pvalue=CollectAsso(pvalue,cutoff,cnv,celltree)
    if (!is.null(pvalue)){
      pvalue$CNA=apply(pvalue,1,function(x,cnv,celltree,cell){
        node=x[2]
        region=x[3]
        case=splitTree(as.character(node),celltree)
        #  control=setdiff(as.character(cell$cell[cell$depth!=0]),case)
        index1=match(case,row.names(cnv))
        #  index2=match(control,row.names(cnv))
        if (mean(cnv[index1[!is.na(index1)],region])>2){
          return("AMP")
        }else if (mean(cnv[index1[!is.na(index1)],region])<2){
          return("DEL")
        }else{
          return(NA)
        }
      },cnv=cnv,celltree=celltree,cell=cell)
      pvalue=pvalue[!is.na(pvalue$CNA),]
      return(pvalue)
    }
  }
}

mergeCNA <- function(node,sigCNA,band.region,arm.band,refer.band){
  refineSigregion <- function(chr,psub,type,refer.band,arm.band,subres){
    if (dim(psub)[1]==1){
      subres=rbind(subres,psub)
    }else if (dim(psub)[1]>1){
      index=match(as.character(psub$region),as.character(refer.band$ID))
      subrefer=refer.band[index,]
      sublength=sum(subrefer$V3-subrefer$V2)
      f=sublength/arm.band$length[arm.band$chr==i&arm.band$band=="p"]
      if (f > 0.5){
        subres=rbind(subres,data.frame(depth=unique(psub$depth),cell=unique(as.character(psub$cell)),region=paste(i,type,sep=":"),monotone.test=mean(psub$monotone.test),trend.test=mean(psub$trend.test),CNA=unique(as.character(psub$CNA))))
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
              return(data.frame(depth=unique(sub1$depth),cell=unique(as.character(sub1$cell)),region=paste(i,x,sep=":"),monotone.test=mean(sub1$monotone.test),trend.test=mean(sub1$trend.test),CNA=unique(as.character(sub1$CNA))))
            }else{
              return(sub1)
            }
          }
        },psub))
        subres=rbind(subres,sub1)
      }
    }
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
      if (length(unique(as.character(psub$CNA)))==1){
        subres=refineSigregion(i,psub,type="p",refer.band,arm.band,subres)
      }else{
        DEL=psub[psub$CNA=="DEL",]
        AMP=psub[psub$CNA=="AMP",]
        subres=refineSigregion(i,DEL,type="p",refer.band,arm.band,subres)
        subres=refineSigregion(i,AMP,type="p",refer.band,arm.band,subres)
      }
    }
    if (dim(qsub)[1]>0){
      if (length(unique(as.character(qsub$CNA)))==1){
        subres=refineSigregion(i,qsub,type="q",refer.band,arm.band,subres)
      }else{
        DEL=qsub[qsub$CNA=="DEL",]
        AMP=qsub[qsub$CNA=="AMP",]
        subres=refineSigregion(i,DEL,type="q",refer.band,arm.band,subres)
        subres=refineSigregion(i,AMP,type="q",refer.band,arm.band,subres)
      }
    }
  }
  return(subres)
}
CombineRegion <- function(node,newsig,refer.band){
  cellsig=c()
  subsig=newsig[newsig$cell==node,]
  indexregion=match(subsig$region,refer.band$ID)
  if (length(indexregion[!is.na(indexregion)])>0){
    subsig1=subsig[!is.na(indexregion),]
    if (dim(subsig1)[1]==1){
      cellsig=rbind(cellsig,subsig1[,1:6])
    }else if (dim(subsig)[1]>1){
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
        if (subsig1$band[k+1]-subsig1$band[k]==1&subsig1$CNA[k+1]==subsig1$CNA[k]&subsig1$arm[k+1]==subsig1$arm[k]){
          k=k+1
          if (k == dim(subsig1)[1]){
            end=k
            region=paste(subsig1$chr[start],":",subsig1$arm[start],subsig1$bandID[start],"-",subsig1$bandID[end],sep="")
            tempsig=data.frame(depth=subsig1$depth[start],cell=subsig1$cell[start],region=region,monotone.test=mean(subsig1$monotone.test[start:end]),trend.test=mean(subsig1$trend.test[start:end]),CNA=subsig1$CNA[start])
            newsig11=rbind(newsig11,tempsig)
          }
        }else{
          end=k
          if (start == end){
            newsig11=rbind(newsig11,subsig1[k,1:6])
          }else{
            region=paste(subsig1$chr[start],":",subsig1$arm[start],subsig1$bandID[start],"-",subsig1$bandID[end],sep="")
            tempsig=data.frame(depth=subsig1$depth[start],cell=subsig1$cell[start],region=region,monotone.test=mean(subsig1$monotone.test[start:end]),trend.test=mean(subsig1$trend.test[start:end]),CNA=subsig1$CNA[start])
            newsig11=rbind(newsig11,tempsig)
          }
          k=k+1
          start=k
          if (start==dim(subsig1)[1]){
            newsig11=rbind(newsig11,subsig1[k,1:6])
          }
        }
      }
      cellsig=rbind(cellsig,newsig11)
    }
  }
  if (length(indexregion[is.na(indexregion)])>0){
    subsig2=subsig[is.na(indexregion),]
    if (dim(subsig2)[1]==1){
      cellsig=rbind(cellsig,subsig2[,1:6])
    }else if (dim(subsig2)[1]>1){
      subregion2=do.call(rbind,sapply(as.character(subsig2$region),strsplit,split=":"))
      index=match(subregion2[,2],c("p","q"))
      if (length(index[!is.na(index)])>0){
        subsig21=subsig2[!is.na(index),]
        cellsig=rbind(cellsig,subsig21[,1:6])
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
          if (subsig22$chr[k+1]==subsig22$chr[k]&subsig22$arm[k+1]==subsig22$arm[k]&subsig22$CNA[k+1]==subsig22$CNA[k]&subsig22$band[k+1]-subsig22$band[k]==1){
            k=k+1
            if (k == dim(subsig22)[1]){
              end=k
              region=paste(subsig22$chr[start],":",subsig22$arm[start],subsig22$band[start],"-",subsig22$band[end],sep="")
              tempsig=data.frame(depth=subsig22$depth[start],cell=subsig22$cell[start],region=region,monotone.test=mean(subsig22$monotone.test[start:end]),trend.test=mean(subsig22$trend.test[start:end]),CNA=subsig22$CNA[start])
              newsig22=rbind(newsig22,tempsig)
            }
          }else{
            end=k
            if (start == end){
              newsig22=rbind(newsig22,subsig22[k,1:6])
            }else{
              region=paste(subsig22$chr[start],":",subsig22$arm[start],subsig22$band[start],"-",subsig22$band[end],sep="")
              tempsig=data.frame(depth=subsig22$depth[start],cell=subsig22$cell[start],region=region,monotone.test=mean(subsig22$monotone.test[start:end]),trend.test=mean(subsig22$trend.test[start:end]),CNA=subsig22$CNA[start])
              newsig22=rbind(newsig22,tempsig)
            }
            k=k+1
            start=k
            if (start==dim(subsig22)[1]){
              newsig22=rbind(newsig22,subsig22[k,1:6])
            }
          }
        }
        cellsig=rbind(cellsig,newsig22)
      }
    }
  }
  return(cellsig)
}
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


ExcludeCNA <- function(cellsig,refer.band,arm.band,band.region,celltree){
  cellsig$index=1
  cellsig$chr=do.call(rbind,strsplit(as.character(cellsig$region),split=":"))[,1]
  cellsig$region1=do.call(rbind,strsplit(as.character(cellsig$region),split=":"))[,2]
  bandID=do.call(rbind,strsplit(as.character(cellsig$region1),split="-"))
  cellsig$start=bandID[,1]
  cellsig$end=apply(bandID,1,function(x){
    if (!substr(x[2],1,1) %in% c("p","q")){
      return(paste(substr(x[1],1,1),x[2],sep=""))
    }else{
      return(x[2])
    }
  })
  index=match(paste(cellsig$chr,":",cellsig$start,sep=""),as.character(refer.band$ID))
  cellsig$startpos[!is.na(index)]=refer.band$V2[index[!is.na(index)]]
  index=match(paste(cellsig$chr,":",cellsig$end,sep=""),as.character(refer.band$ID))
  cellsig$endpos[!is.na(index)]=refer.band$V3[index[!is.na(index)]]
  index=match(paste(cellsig$chr,":",cellsig$start,sep=""),paste(arm.band$chr,arm.band$band,sep=":"))
  cellsig$startpos[!is.na(index)]=arm.band$start[index[!is.na(index)]]
  index=match(paste(cellsig$chr,":",cellsig$end,sep=""),paste(arm.band$chr,arm.band$band,sep=":"))
  cellsig$endpos[!is.na(index)]=arm.band$end[index[!is.na(index)]]
  index=match(paste(cellsig$chr,":",cellsig$start,sep=""),paste(band.region$chr,band.region$band,sep=":"))
  cellsig$startpos[!is.na(index)]=band.region$start[index[!is.na(index)]]
  index=match(paste(cellsig$chr,":",cellsig$end,sep=""),paste(band.region$chr,band.region$band,sep=":"))
  cellsig$endpos[!is.na(index)]=band.region$end[index[!is.na(index)]]
  maxD=max(as.numeric(as.character(cellsig$depth)))
  for (i in 1:dim(cellsig)[1]){
    if (as.numeric(as.character(cellsig$depth[i]))<maxD){
      child=splitTree(as.character(cellsig$cell[i]),celltree)
      child=setdiff(child,as.character(cellsig$cell[i]))
      start=cellsig$startpos[i]
      end=cellsig$endpos[i]
      chr=cellsig$chr[i]
      index=match(child,as.character(cellsig$cell))
      subcell=cellsig[index[!is.na(index)],]
      subcell=subcell[subcell$chr==chr,]
      if (dim(subcell)[1]>0){
        deltStart=subcell$startpos-start
        deltEnd=subcell$endpos-end
        k=which(deltStart>=0&deltEnd<=0)
        if (length(k)>0){
          cellsig$index[as.character(cellsig$cell)==as.character(subcell$cell)[k]&as.character(cellsig$region)==as.character(subcell$region)[k]]=0
        }
      }
    }
  }
  return(cellsig$index)
}


LSAdefine <- function(sigCNA,datapath){
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
  newsig=do.call(rbind,lapply(unique(as.character(sigCNA$cell)),mergeCNA,sigCNA,band.region,arm.band,refer.band))
  newsig$event=paste(as.character(newsig$region),as.character(newsig$CNA),sep="_")
  event=table(newsig$event)
  if (length(event[event>1])>0){
    event=event[event>1]
    for (j in names(event)){
      subsig=newsig[newsig$event==j,]
      overlap=unlist(lapply(subsig$cell,function(x,celltree,subsig){
        return(intersect(setdiff(splitTree(x,celltree),x),subsig$cell))
      },celltree,subsig))
      if (length(overlap)>0){
        index=match(paste(overlap,j,sep=':'),paste(as.character(newsig$cell),as.character(newsig$event),sep=":"))
        newsig=newsig[-index,]
      }
    }
  }
  cell=unique(as.character(newsig$cell))
  cellsig=do.call(rbind,lapply(cell,CombineRegion,newsig,refer.band))
  cellsig$index=ExcludeCNA(cellsig,refer.band,arm.band,band.region,celltree)
  cellsig1=cellsig[cellsig$index==1,]
  return(cellsig1)
}



CNAconnect <- function(sigCNA,celltree){
  node=unique(sigCNA[,1:2])
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
        node1=as.character(node[!is.na(index),][1,2])
        dist=as.numeric(as.character(node[i,1]))-as.numeric(as.character(node[!is.na(index),][1,1]))
      }
      CNAnetwork=rbind(CNAnetwork,c(node1,as.character(node$cell[i]),dist))
    }
  }
  colnames(CNAnetwork)=c("node1","node2","dist")
  return(CNAnetwork)
}


geneCNAfunction<-function(k,pathwaygene,ancestorCNV,generegion){
  gene=pathwaygene[k,]
  chr=as.character(gene[1,1])
  start=gene[1,2]
  end=gene[1,3]
  gene=gene[1,4]
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
celltree=read.csv(treeName,sep="\t")
nodes=data.frame(id=union(as.character(celltree[,1]),as.character(celltree[,2])),size=5)
nodes$color="lightblue"
nodes$color[nodes$id==setdiff(as.character(celltree[,1]),as.character(celltree[,2]))]="black"
net <- graph_from_data_frame(d=celltree, vertices=nodes, directed=T)
pdf(file=paste(outpath,"/singlecell.tree.pdf",sep=""),width = 5,height = 5,useDingbats = F)
plot(net, vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label=NA)
dev.off()
###lineage specific segmental CNAs
data = read.csv(inputfile,sep="\t",header = TRUE)
if (datatype=="D"){
  region=data[,1:2]
  region$chrom=paste("chr",region$chrom,sep="")
  region$end=region$chrompos
}else if (datatype=="R"){

  reference=paste(datapath,"/gencode_v19_gene_pos.txt",sep="")
  geneInfo=read.csv(reference,sep="\t",header=F)
  index=match(row.names(data),as.character(geneInfo[,1]))
  newdata=cbind(geneInfo[index[!is.na(index)],2:3],data[!is.na(index),])
  rownames(newdata)=rownames(data)[!is.na(index)]
  #ll=nchar(as.character(newdata[,1]))
  #chromo=data.frame(chr=as.character(newdata[,1]),ll=ll)
  #chro=apply(chromo,1,function(x){
  #  return(substr(x[1],start=4,stop=x[2]))
  #})
  #chro[chro=="chrX"]="chr23"
  #newdata[,1]=chro
  newdata=newdata[as.character(newdata[,1])!="chrM"&as.character(newdata[,1])!="chrY",]
  newdata[,1]=as.character(newdata[,1])
  newdata[newdata[,1]=="chrX",1]="chr23"
  region=newdata[,1:2]
  #region$chrom=paste("chr",region$chrom,sep="")
  region[,3]=region[,2]
  colnames(region)=c("chrom","chrompos","end")
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

cell=union(as.character(celltree[,1]),as.character(celltree[,2]))
cell=data.frame(cell=cell)
cell$depth=sapply(as.character(cell$cell),depthFunction,cellTree=celltree)
cell$subtreesize=sapply(as.character(cell$cell),subtreeSize,cellTree=celltree)
sigCNA=RefineCNA(cnv=newCNV,cell,celltree,cutoff = 0.01)
if (!is.null(dim(sigCNA))){
  sigCNA=sigCNA[order(sigCNA$depth),]
  sigCNA=LSAdefine(sigCNA,datapath)
  sigCNA=sigCNA[,1:(dim(sigCNA)[2]-1)]
  write.table(sigCNA,paste(outpath,"/segmental.LSA.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}else{
  print.noquote("No segmental LSA is identified!")
}

####lineage specific cancer associated genes
pathwaygene=read.csv(paste(datapath,"/pathwaygene.txt",sep=""),sep="\t")
if (datatype=="D"){
  cnv=read.csv(CNVfile,sep="\t")
  oncogenicCNV=do.call(cbind,lapply(1:dim(pathwaygene)[1],geneCNAfunction,pathwaygene=pathwaygene,ancestorCNV=cnv,generegion=region))
  colnames(oncogenicCNV)=as.character(pathwaygene$name)
}else if (datatype == "R"){
  data=data[,3:dim(data)[2]]
  index=match(as.character(pathwaygene$name),rownames(data))
  oncogenicCNV=t(data[index[!is.na(index)],])
  oncogenicCNV=round(oncogenicCNV*2)
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
siggene=RefineCNA(cnv=oncogenicCNV,cell,celltree,cutoff = 0.01)
if (!is.null(dim(siggene))){
  write.table(siggene,paste(outpath,"/gene.LSA.txt",sep=""),col.names = T,row.names=T,sep="\t",quote=F)
}else{
  print.noquote("No LSA is identified at gene level!")
}


####plot LSA Tree
allsig=c()
if (!is.null(dim(sigCNA))){
  allsig=sigCNA
  if (!is.null(dim(siggene))){
    allsig=rbind(allsig,siggene)
  }
}else{
  if (!is.null(dim(siggene))){
    allsig=siggene
  }
}
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
      pvalue=apply(CNA[,4:5],1,min)
      CNA=CNA[order(pvalue),]
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
  links=data.frame(from=LSAnetwork[,1],to=LSAnetwork[,2],weight=as.numeric(LSAnetwork[,3]))
  pdf(file=paste(outpath,"/LSA.tree.pdf",sep=""),width = 5,height = 5,useDingbats = F)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
  plot(net, layout=layout_as_tree,vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label.cex=0.5,vertex.label=nodes$annotation)
  dev.off()
}
