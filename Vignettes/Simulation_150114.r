

rm(list=ls(all=T))

library(MASS)
library(foreach);library(doSNOW)
library(Rcpp)

# library(tCoxnet)
source("CV_cutoff.R")
source("loCoxNet11.R")
sourceCpp("loCoxNet11.cpp")
# 

########################
#####  Simulation  #####
########################

#####  Data parameters  #####

fb1=function(x){2/(exp(x^3))}
fb2=function(x){4*x*(1-x)}
fb3=function(x){(sin(pi*x)+cos(pi*x))}

fb=function(x,index=c(1,2,3)){
  out=NULL
  for(i in 1:length(index)){
    ii=abs(index[i]);si=sign(index[i])
    out=c(out,si*switch(ii,"1"=fb1(x),"2"=fb2(x),"3"=fb3(x)))}
  return(out)
}

N=500;r=0.2
ngroup=c(3,20);nb1=rep(5,ngroup[1]);nb2=rep(5,ngroup[2])
p1=sum(nb1);p2=sum(nb2);p=p1+p2
model_index=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)
#model_index=c(1,1,1,-1,-1,2,2,2,-2,-2,3,3,3,-3,-3)

a1=5;a2=2

### True adjacent matrix
Omega0=matrix(0,ncol=p,nrow=p);ij=1
for(i in 1:ngroup[1]){
  R=matrix(r,ncol=nb1[i],nrow=nb1[i]);diag(R)=1
  Omega0[ij:(ij+nb1[i]-1),ij:(ij+nb1[i]-1)]=R
  ij=ij+nb1[i]
}
for(i in 1:ngroup[2]){
  R=matrix(r,ncol=nb2[i],nrow=nb2[i]);diag(R)=1
  Omega0[ij:(ij+nb2[i]-1),ij:(ij+nb2[i]-1)]=R
  ij=ij+nb2[i]
}
Adjm0=Omega0;diag(Adjm0)=0;Adjm0[Adjm0>0]=1


#####  Simulation  #####

###  algorithm parameters  ###
w0=seq(0.1,0.9,by=0.1);nw0=length(w0)
hs=c(0.1,0.2,0.3,0.4);nh=length(hs)
b0=do.call("cbind",lapply(w0,function(x){c(fb(x,model_index),rep(0,p2))}))

di=-0.01

alphas=seq(0.1,0.9,by=0.1);nalpha=length(alphas)
lambda=NULL;nlambda=20
nfolds=10;alinear=TRUE;alambda=NULL;nalambda=10

ncluster=4
setDefaultClusterOptions(type="SOCK")
getClusterOption("type")
getDoParWorkers()

###  simulation parameters  ###
set.seed(1213)

Beta=list();Betat=list();Betase=list()
select.para=NULL;cuts=NULL
Iterations=100;ite=1

###  Run  ###
for(ite in 1:Iterations){
  cat(ite);write.table(ite,"Current.txt",row.names=F,col.names=F)
  betai=list();betat=list();betase=list()
  
  #####  Data generation  #####
  
  ###  Step 1 - w,beta,x  ###
  w=runif(N,0,1)
  b1=do.call("rbind",lapply(w,function(x){fb(x,model_index)}))
  
  ## x: relevant
  x=NULL
  for(i in 1:ngroup[1]){
    R=matrix(r,ncol=nb1[i],nrow=nb1[i]);diag(R)=1
    x=cbind(x,mvrnorm(N,mu=rep(0,nb1[i]),Sigma=R))
  }
  ## x: irrelevant
  for(i in 1:ngroup[2]){
    R=matrix(r,ncol=nb2[i],nrow=nb2[i]);diag(R)=1
    x=cbind(x,mvrnorm(N,mu=rep(0,nb2[i]),Sigma=R))
  }
  
  #x=apply(x,2,function(tx){tx-mean(tx)})
  x=apply(x,2,function(tx,nr){sdn=sd(tx)*sqrt((nr-1)/nr);return((tx-mean(tx))/sdn)},nr=nrow(x))
  
  Omega=cor(x);Omega=abs(Omega)
  
  ###  Step 2 - time and censoring  ###
  u=runif(N,0,1)
  xb=apply(x[,1:p1]*b1,1,sum)
  ty=(-log(u)/a2/exp(xb))^(1/a1)
  tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
  #tcens=as.numeric(runif(N,2,15)<=ty)
  y=cbind(time=ty,status=1-tcens)
  
  
  #####  Initial  #####
  
  ###  Split data at the outset  ###
  foldid=split.foldidw(w,nrow(x),nfolds)
  
  
  #####  Estimation  #####
  
  save(list=ls(all=T),file="tem.rdata")
  
  ###  Model 1 - AdalocNet+Omega0  ###
  
  ###  Intial values for adaptive estimates  ###
  modeli=0
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  aini=list()
  for (ih in 1:nh) {
    aini[[ih]]=locoxini(x,y,w,w0,hs[ih],alambda,nalambda,rlambda=NULL,alinear,1,foldid,N,p)
  }
  # stopCluster(cll)
  
  
  modeli=1;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outj=list()
  for(ih in 1:nh){
    ia=1;outi=list()
    repeat{
      outi[[ia]]=ttCoxnet(x,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Network",Omega=Omega0,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini[[ih]],nfolds=nfolds,foldid=foldid,keep.beta=T)
      if(ia>1){if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm))break}
      if(ia==nalpha)break
      ia=ia+1
    }
    if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm)){
      cv.fiti=outi[[ia-1]];iai=ia-1
    }else{
      cv.fiti=outi[[ia]];iai=ia
    }
    outj[[ih]]=(list(fit=cv.fiti,ia=iai))
  }
  # stopCluster(cll)
  
  ### Max
  index.max=sapply(outj,function(x){which.max(x$fit$fit$cvm)})
  
  ## Max-Beta
  tembetai=list()
  for(ih in 1:nh){tembetai[[ih]]=outj[[ih]]$fit$Beta[[index.max[ih]]]}
  betai[[modeli]]=tembetai
  
  ## Max-para
  temparai=NULL
  for(ih in 1:nh){
    temparai=data.frame(alpha=alphas[outj[[ih]]$ia],lambda=outj[[ih]]$fit$fit$lambda[index.max[ih]],h=hs[ih],cvh0=outj[[ih]]$fit$cvh$cvm[1],cvh1=outj[[ih]]$fit$cvh$cvm[2],model=modeli,index=ite)
    select.para=rbind(select.para,temparai)
  }
  
  ## Max-h
  dcvhi=numeric(nh)
  for(ih in 1:nh){dcvhi[ih]=outj[[ih]]$fit$cvh$cvm[1]-outj[[ih]]$fit$cvh$cvm[2]}
  tem=which((dcvhi)>di)
  if(length(tem)>0){indexh=tem[1]}else{indexh=nh}
  # indexh;mse.bandwidth2(outj,index.max,b0)
  
  cv.fiti0=outj[[indexh]]$fit
  indexi=which.max(cv.fiti0$fit$cvm)
  indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))
  indexij=which(cv.fiti0$fit$index=="se")
  
  ### Trimming
  betati=list()
  ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
  absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)
  absmj=apply(abs(cv.fiti0$Beta[[indexj]]),1,mean)
  
  #Trim1
  cutoffj=seq(0,max(absmi)-0.001,length=20)  
  cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,adaptive=T,aini=aini[[indexh]],nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,adaptive=T,aini=aini[[indexh]],nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi1=cuti$cutoff[ddi0]
  
  ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0
  
  #Trim2
  cutoffj=seq(0,max(absmj)-0.001,length=20)  
  cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffj,adaptive=T,aini=aini[[indexh]],nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffi,adaptive=T,aini=aini[[indexh]],nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi2=cuti$cutoff[ddi0]
  
  ib=3;betati[[ib]]=cv.fiti0$Beta[[indexj]];betati[[ib]][absmj<=cutoffi2,]=0
  
  ib=4;betati[[ib]]=cv.fiti0$Beta[[indexij]]
  
  betat[[modeli]]=betati
  cutsi=data.frame(cut1=cutoffi1,cut2=cutoffi2,model=modeli,index=ite)
  cuts=rbind(cuts,cutsi)
  
  
  save(list=ls(all=T),file="tem.rdata")
  
  
  
  ###  Model 2 - locNet+Omega0  ###
  modeli=2;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outj=list()
  for(ih in 1:nh) {
    ia=1;outi=list()
    repeat{
      outi[[ia]]=ttCoxnet(x,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Network",Omega=Omega0,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
      if(ia>1){if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm))break}
      if(ia==nalpha)break
      ia=ia+1
    }
    if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm)){
      cv.fiti=outi[[ia-1]];iai=ia-1
    }else{
      cv.fiti=outi[[ia]];iai=ia
    }
    outj[[ih]]=(list(fit=cv.fiti,ia=iai))
  }
  # stopCluster(cll)
  
  ### Max
  index.max=sapply(outj,function(x){which.max(x$fit$fit$cvm)})
  
  ## Max-Beta
  tembetai=list()
  for(ih in 1:nh){tembetai[[ih]]=outj[[ih]]$fit$Beta[[index.max[ih]]]}
  betai[[modeli]]=tembetai
  
  ## Max-para
  temparai=NULL
  for(ih in 1:nh){
    temparai=data.frame(alpha=alphas[outj[[ih]]$ia],lambda=outj[[ih]]$fit$fit$lambda[index.max[ih]],h=hs[ih],cvh0=outj[[ih]]$fit$cvh$cvm[1],cvh1=outj[[ih]]$fit$cvh$cvm[2],model=modeli,index=ite)
    select.para=rbind(select.para,temparai)
  }
  
  ## Max-h
  dcvhi=numeric(nh)
  for(ih in 1:nh){dcvhi[ih]=outj[[ih]]$fit$cvh$cvm[1]-outj[[ih]]$fit$cvh$cvm[2]}
  tem=which((dcvhi)>di)
  if(length(tem)>0){indexh=tem[1]}else{indexh=nh}
  # indexh;mse.bandwidth2(outj,index.max,b0)
  
  cv.fiti0=outj[[indexh]]$fit
  indexi=which.max(cv.fiti0$fit$cvm)
  indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))
  indexij=which(cv.fiti0$fit$index=="se")
  
  ### Trimming
  betati=list()
  ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
  absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)
  absmj=apply(abs(cv.fiti0$Beta[[indexj]]),1,mean)
  
  #Trim1
  cutoffj=seq(0,max(absmi)-0.001,length=20)  
  cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi1=cuti$cutoff[ddi0]
  
  ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0
  
  #Trim2
  cutoffj=seq(0,max(absmj)-0.001,length=20)  
  cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],Omega=Omega0,alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi2=cuti$cutoff[ddi0]
  
  ib=3;betati[[ib]]=cv.fiti0$Beta[[indexj]];betati[[ib]][absmj<=cutoffi2,]=0
  
  ib=4;betati[[ib]]=cv.fiti0$Beta[[indexij]]
  
  betat[[modeli]]=betati
  cutsi=data.frame(cut1=cutoffi1,cut2=cutoffi2,model=modeli,index=ite)
  cuts=rbind(cuts,cutsi)
  
  
  save(list=ls(all=T),file="tem.rdata")
  
  
  
  ###  Model 3 - locEnet  ###
  modeli=3;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outj=list()
  for(ih in 1:nh) {
    ia=1;outi=list()
    repeat{
      outi[[ia]]=ttCoxnet(x,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Enet",alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
      if(ia>1){if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm))break}
      if(ia==nalpha)break
      ia=ia+1
    }
    if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm)){
      cv.fiti=outi[[ia-1]];iai=ia-1
    }else{
      cv.fiti=outi[[ia]];iai=ia
    }
    outj[[ih]]=(list(fit=cv.fiti,ia=iai))
  }
  # stopCluster(cll)
  
  ### Max
  index.max=sapply(outj,function(x){which.max(x$fit$fit$cvm)})
  
  ## Max-Beta
  tembetai=list()
  for(ih in 1:nh){tembetai[[ih]]=outj[[ih]]$fit$Beta[[index.max[ih]]]}
  betai[[modeli]]=tembetai
  
  ## Max-para
  temparai=NULL
  for(ih in 1:nh){
    temparai=data.frame(alpha=alphas[outj[[ih]]$ia],lambda=outj[[ih]]$fit$fit$lambda[index.max[ih]],h=hs[ih],cvh0=outj[[ih]]$fit$cvh$cvm[1],cvh1=outj[[ih]]$fit$cvh$cvm[2],model=modeli,index=ite)
    select.para=rbind(select.para,temparai)
  }
  
  ## Max-h
  dcvhi=numeric(nh)
  for(ih in 1:nh){dcvhi[ih]=outj[[ih]]$fit$cvh$cvm[1]-outj[[ih]]$fit$cvh$cvm[2]}
  tem=which((dcvhi)>di)
  if(length(tem)>0){indexh=tem[1]}else{indexh=nh}
  # indexh;mse.bandwidth2(outj,index.max,b0)
  
  cv.fiti0=outj[[indexh]]$fit
  indexi=which.max(cv.fiti0$fit$cvm)
  indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))
  indexij=which(cv.fiti0$fit$index=="se")
  
  ### Trimming
  betati=list()
  ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
  absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)
  absmj=apply(abs(cv.fiti0$Beta[[indexj]]),1,mean)
  
  #Trim1
  cutoffj=seq(0,max(absmi)-0.001,length=20)  
  cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi1=cuti$cutoff[ddi0]
  
  ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0
  
  #Trim2
  cutoffj=seq(0,max(absmj)-0.001,length=20)  
  cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],alpha=alphas[outj[[indexh]]$ia],lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi2=cuti$cutoff[ddi0]
  
  ib=3;betati[[ib]]=cv.fiti0$Beta[[indexj]];betati[[ib]][absmj<=cutoffi2,]=0
  
  ib=4;betati[[ib]]=cv.fiti0$Beta[[indexij]]
  
  betat[[modeli]]=betati
  cutsi=data.frame(cut1=cutoffi1,cut2=cutoffi2,model=modeli,index=ite)
  cuts=rbind(cuts,cutsi)
  
  
  save(list=ls(all=T),file="tem.rdata")
  
  
  
  ###  Model 4 - locL1  ###
  modeli=4;alpha=1.0;lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  cll=makeCluster(ncluster);registerDoSNOW(cll)
  outj=list()
  for(ih in 1:nh) {
    cv.fiti=ttCoxnet(x,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Enet",alpha=alpha,lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
    outj[[ih]]=(cv.fiti)
  }
  stopCluster(cll)
  
  ### Max
  index.max=sapply(outj,function(x){which.max(x$fit$cvm)})
  
  ## Max-Beta
  tembetai=list()
  for(ih in 1:nh){tembetai[[ih]]=outj[[ih]]$Beta[[index.max[ih]]]}
  betai[[modeli]]=tembetai
  
  ## Max-para
  temparai=NULL
  for(ih in 1:nh){
    temparai=data.frame(alpha=alpha,lambda=outj[[ih]]$fit$lambda[index.max[ih]],h=hs[ih],cvh0=outj[[ih]]$cvh$cvm[1],cvh1=outj[[ih]]$cvh$cvm[2],model=modeli,index=ite)
    select.para=rbind(select.para,temparai)
  }
  
  ## Max-h
  dcvhi=numeric(nh)
  for(ih in 1:nh){dcvhi[ih]=outj[[ih]]$cvh$cvm[1]-outj[[ih]]$cvh$cvm[2]}
  tem=which((dcvhi)>di)
  if(length(tem)>0){indexh=tem[1]}else{indexh=nh}
  # indexh
  
  cv.fiti0=outj[[indexh]]
  indexi=which.max(cv.fiti0$fit$cvm)
  indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))
  indexij=which(cv.fiti0$fit$index=="se")
  
  ### Trimming
  betati=list()
  ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
  absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)
  absmj=apply(abs(cv.fiti0$Beta[[indexj]]),1,mean)
  
  #Trim1
  cutoffj=seq(0,max(absmi)-0.001,length=20)  
  cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexi]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi1=cuti$cutoff[ddi0]
  
  ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0
  
  #Trim2
  cutoffj=seq(0,max(absmj)-0.001,length=20)  
  cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cutj$cutoff)
  # cbind(cutj,mseii)
  # plot(cutj$cutoff,cutj$cvm,type="b")
  ddi0=ddj0=5
  repeat{
    dcvm=max(cutj$cvm[ddi0:ddj0])-min(cutj$cvm[ddi0:ddj0])
    ddij=which.max(cutj$cvm[ddi0:ddj0])+ddi0-1
    if(dcvm>cutj$cvse[ddij])break
    if(ddj0>15){ddj0=10;break}
    ddj0=ddj0+1
  }
  cutoffi=seq(0,cutj$cutoff[ddj0],length=20)
  cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexj]],x,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexj],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
  # mseii=mse.cutoff(cv.fiti0$Beta[[indexj]],b0,cuti$cutoff)
  # cbind(cuti,mseii)
  # plot(cuti$cutoff,cuti$cvm,type="b")
  ddi0=15
  repeat{
    if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
    if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
    ddi0=ddi0-1
  }
  cutoffi2=cuti$cutoff[ddi0]
  
  ib=3;betati[[ib]]=cv.fiti0$Beta[[indexj]];betati[[ib]][absmj<=cutoffi2,]=0
  
  ib=4;betati[[ib]]=cv.fiti0$Beta[[indexij]]
  
  betat[[modeli]]=betati
  cutsi=data.frame(cut1=cutoffi1,cut2=cutoffi2,model=modeli,index=ite)
  cuts=rbind(cuts,cutsi)
  
  
  save(list=ls(all=T),file="tem.rdata")
  
  
  
  ###  Method 5 - AdaNet+Omega0  ###
  
  ### Intial values
  modeli=10
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  outi=ttCoxnet(x,y,penalty="Enet",alpha=0.0,lambda=0.0) # MLE
  aini2=list(wbeta=1/abs(outi$Beta),sgn=sign(outi$Beta))
  
  ###
  modeli=5;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outi=list()
  for(ia in 1:nalpha) {
    fiti=ttCoxnet(x,y,penalty="Network",Omega=Omega0,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini2,nfolds=nfolds,foldid=foldid,keep.beta=T)
    outi[[ia]]=(fiti)
  }
  # stopCluster(cll)
  
  ### Max
  indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
  indexi=which.max(outi[[indexa]]$fit$cvm)
  indexij=which(outi[[indexa]]$fit$index=="se")
  
  temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli,index=ite)
  select.para=rbind(select.para,temparai)
  
  betai[[modeli]]=outi[[indexa]]$Beta[,indexi]
  betase[[modeli]]=outi[[indexa]]$Beta[,indexij]
  
  
  ###  Method 6 - Net+Omega0  ###
  modeli=6;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outi=list()
  for(ia in 1:nalpha) {
    fiti=ttCoxnet(x,y,penalty="Network",Omega=Omega0,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
    outi[[ih]]=(fiti)
  }
  # stopCluster(cll)
  
  ### Max
  indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
  indexi=which.max(outi[[indexa]]$fit$cvm)
  indexij=which(outi[[indexa]]$fit$index=="se")
  
  temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli,index=ite)
  select.para=rbind(select.para,temparai)
  
  betai[[modeli]]=outi[[indexa]]$Beta[,indexi]
  betase[[modeli]]=outi[[indexa]]$Beta[,indexij]
  
  
  ###  Model 7 - Enet  ###
  modeli=7;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  # cll=makeCluster(ncluster);registerDoSNOW(cll)
  outi=list()
  for(ia in 1:nalpha) {
    fiti=ttCoxnet(x,y,penalty="Enet",alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
    outi[[ia]]=(fiti)
  }
  # stopCluster(cll)
  
  ### Max
  indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
  indexi=which.max(outi[[indexa]]$fit$cvm)
  indexij=which(outi[[indexa]]$fit$index=="se")
  
  temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli,index=ite)
  select.para=rbind(select.para,temparai)
  
  betai[[modeli]]=outi[[indexa]]$Beta[,indexi]
  betase[[modeli]]=outi[[indexa]]$Beta[,indexij]
  
  
  ###  Method 8 - L1  ###
  modeli=8;alpha=1.0;lambda=NULL
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  fiti=ttCoxnet(x,y,penalty="Enet",alpha=alpha,lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
  
  indexi=which.max(fiti$fit$cvm)
  indexij=which(fiti$fit$index=="se")
  
  temparai=data.frame(alpha=alpha,lambda=fiti$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli,index=ite)
  select.para=rbind(select.para,temparai)
  
  betai[[modeli]]=fiti$Beta[,indexi]
  betase[[modeli]]=fiti$Beta[,indexij]
  
  
  ###  Method 9 - None  ###
  modeli=9;alpha=1.0;lambda=0.0
  write.table(modeli,"Model.txt",row.names=F,col.names=F)
  
  fiti=ttCoxnet(x,y,penalty="Enet",alpha=alpha,lambda=lambda)
  
  betai[[modeli]]=fiti$Beta
  
  
  ###  record  ###
  Beta[[ite]]=betai;Betat[[ite]]=betat;Betase[[ite]]=betase
  
  save(Beta,Betat,Betase,select.para,cuts,file="Simulation_results.rdata")
}

save(Beta,Betat,Betase,select.para,cuts,file="Simulation_results.rdata")



b0=do.call("cbind",lapply(w0,function(x){c(fb(x,model_index),rep(0,p2))}))
result=list();result2=list()
for(ij in 1:2){
  i=1;mse=NULL;nzeros=NULL
  for(i in 1:length(Betat[[ij]])){
    tembeta=Betat[[ij]][[i]];msei=numeric(4);nzeroi=numeric(4)
    for(j in 1:4){
      msei[j]=sum(apply((tembeta[[j]]-b0)^2,1,mean))
      abmi=apply(abs(tembeta[[j]]),1,mean)
      nzeroi[j]=sum(abmi!=0)
    }
    mse=rbind(mse,msei);nzeros=rbind(nzeros,nzeroi)
  }
  mse;nzeros
  result[[ij]]=mse;result2[[ij]]=nzeros
}


mse=NULL
for(i in 5:9){mse=c(mse,sum(apply((betai[[i]]-b0)^2,1,mean)))}
mse


