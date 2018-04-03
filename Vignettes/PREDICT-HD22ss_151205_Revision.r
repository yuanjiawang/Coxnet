


################################
#####  Real Data analysis  #####
################################

rm(list=ls(all=T))

load("Data.rdata")

library(car)
source("CV_cutoff.R")
source("loCoxNet11.R")
sourceCpp("loCoxNet11.cpp")
library(glmnet)
# library(foreach);library(doSNOW)


di=-0.01

alphas=seq(0.1,0.9,by=0.1);nalpha=length(alphas)
lambda=NULL;nlambda=50
nfolds=10;alinear=TRUE;alambda=NULL;nalambda=10

foldid=split.foldidw(w,nrow(x),nfolds)

ncluster=3
setDefaultClusterOptions(type="SOCK")
getClusterOption("type")
getDoParWorkers()



#####################
#####  AlocNet  #####
#####################
w0=seq(0.50,1.00,by=0.01);nw0=length(w0)

# h=0.15;ia=1
# aini=locoxini(x=xstd,y,w,w0,h,alambda,nalambda,rlambda=NULL,alinear=F)
# ia=1;outi=list();alphas=seq(0.1,0.9,by=0.1);nalpha=length(alphas)
# repeat{
#   outi[[ia]]=ttCoxnet(x=xstd,y=y,w=w,w0=w0,h=h,hnext=c(h,h+0.1),penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid,keep.beta=T)
#   if(ia>1){if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm))break}
#   if(ia==nalpha)break
#   ia=ia+1
# }
# if(max(outi[[ia]]$fit$cvm)<max(outi[[ia-1]]$fit$cvm)){
#   cv.fiti=outi[[ia-1]];iai=ia-1
# }else{
#   cv.fiti=outi[[ia]];iai=ia
# }
# cv.fiti$cvh;iai

# h=0.2;ia=3;
# h=0.175;ia=3;
# h=0.15;ia=1;

h=0.125;ia=1

modeli=1;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
aini=locoxini(x=xstd,y,w,w0,h,alambda,nalambda,rlambda=NULL,alinear=F)
# outj=ttCoxnet(x=xstd,y=y,w=w,w0=w0,h=h,hnext=c(h,h+0.1),penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid,keep.beta=T)
# outj$cvh

sum(duplicated(y[y[,"status"]==1,"time"]))

outj=Coxnet(x=xstd,y=y,w=w,w0=w0,h=h,hnext=c(h,h+0.1),penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid,keep.beta=T,thresh=1e-5)

lambdaj0=outj$lambda.max

### Trimming
cv.fiti0=outj
indexi=which.max(cv.fiti0$fit$cvm)
indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))

betati=list()
ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)

#Trim1
cutoffj=seq(0,max(absmi)-0.001,length=20)
cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=h,Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid)
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
cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=h,Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid)
# plot(cuti$cutoff,cuti$cvm,type="b")
ddi0=15
repeat{
  if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
  if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
  ddi0=ddi0-1
}
cutoffi1=cuti$cutoff[ddi0]
ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0

Betai=betati[[2]]

nzeroj0=cuti$nzero[ddi0]


### output
out=Betai
out=cbind(out,apply(abs(out),1,mean))
out=cbind(out,apply(out,1,mean))
out=cbind(out,rank(exp(-out[,nw0+1])))
out[out==0]=NA
colnames(out)=c(w0,"M","A","R")

write.csv(out,"R_tmp.csv",na="")

#colnames(Betai)=w0;rownames(Betai)=c(namex,name_image)
#write.csv(Betai,"Trace50.csv")

table(out[,"R"])


###################################################
###  Bootstrap - fix all the tuning parameters  ###
###################################################

set.seed(123)
BetaS=list()

h=0.125;ia=1
modeli=1;alphas=seq(0.1,0.9,by=0.1);B=100
for (ib in 1:B) {
  print(ib)
  
  ## re-sampling
  indexi1=sample(which(y[,"status"]==1), replace=T)
  indexi0=sample(which(y[,"status"]==0), replace=T)
  indexi=c(indexi1,indexi0)
  
  xstdi=xstd[indexi,]; yi=y[indexi,]; wi=w[indexi]
  
  # sum(duplicated(yi[yi[,"status"]==1,"time"]))
  
  ##
  #outj2=ttCoxnet(x=xstdi,y=yi,w=wi,w0=w0,h=h,hnext=NULL,penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambdaj0,nlambda=nlambda,adaptive=T,aini=aini,nfolds=1,foldid=foldid,keep.beta=T)
  
  outj=Coxnet(x=xstdi,y=yi,w=wi,w0=w0,h=h,hnext=NULL,penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambdaj0,nlambda=nlambda,adaptive=T,aini=aini,nfolds=1,foldid=foldid,keep.beta=T,thresh=1e-5)
  
  # range(outj2$Beta[[1]]-outj$Beta[[1]])
  #outj$cvh; outj2$cvh
  #cbind(outj$fit, outj2$fit)
  
  ### Trimming
  cv.fiti0=outj
  absmi=apply(abs(cv.fiti0$Beta[[1]]),1,mean)
  
  betati=cv.fiti0$Beta[[1]]#;betati[absmi<sort(absmi,decreasing=T)[nzeroj0],]=0
  BetaS[[ib]]=betati
}

length(BetaS)

# BetaBooT=BetaS # set.seed(123)
# save(BetaBooT,file="BetaBoot1.tie.rdata")
# load("BetaBoot1.tie.rdata")

# BetaBooT2=BetaS # set.seed(123)
# save(BetaBooT2,file="BetaBoot2.tie.rdata")
# load("BetaBoot2.tie.rdata")



####################
###  Plot trace  ###
####################

out=Betai
out=cbind(out,apply(abs(out),1,mean))
out=cbind(out,rank(exp(-out[,nw0+1])))
colnames(out)=c(w0,"M","R")
rownames(out)=c(namex,name_image)

rownames(out)[out[,"R"]<=14]


## Bootstrap results - fixed tuning parameters
# load("BetaBoot1.tie.rdata")
# BetaBoot=BetaBooT

# save(list=ls(),file="BetaBoot2.tie.all.rdata")

load("BetaBoot2.tie.rdata")
BetaBoot=BetaBooT2

###  non-image  ###
namei=c("TMS","SDMT")
coli=c("black","red2")

indexi=which(rownames(out) %in% c("TOTAL_MOTOR_SCORE", "SDMT"))
betai=list(); i=1
for (i in 1:length(indexi)) {
  betai[[i]]=sapply(BetaBoot,function(x){x[indexi[i],]})
}

# 95% CI
CIup=NULL; CIlow=NULL
for (i in 1:length(indexi)) {
  CIup=cbind(CIup,apply(betai[[i]],1,quantile,prob=0.975))
  CIlow=cbind(CIlow,apply(betai[[i]],1,quantile,prob=0.025))
}

#pdf("Beta_nonimage2.pdf")
#tiff("Figure/Beta_nonimage3.tif",width=6,height=6,units='in',res=300)
#tiff(paste("Beta_nonimage",h,".tif",sep=""),width=6,height=6,units='in',res=300)

cairo_ps("Figure_revision/Beta_nonimage2.eps")

plot(w0,out["TOTAL_MOTOR_SCORE",1:nw0],ylim=c(-2,2),type="l",ylab="Beta",xlab=expression("CAP"[s]),main="Clinical markers",cex.main=1.7,cex.lab=1.3)
abline(h=0,col="grey")

points(w0,out["TOTAL_MOTOR_SCORE",1:nw0],type="l",col=coli[1],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,1],rev(CIup[,1])),col="#D4D4D488",border=NA)

points(w0,out["SDMT",1:nw0],type="l",col=coli[2],lty=2,lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,2],rev(CIup[,2])),col="#D4D4D488",border=NA)

abline(v=0.67,lty="dashed",col="gray30",lwd=1.5);abline(v=0.85,lty="dashed",col="gray30",lwd=1.5)
legend("topleft",legend=namei,col=coli,text.col=coli,lty=c(1,2),lwd=2,cex=1.3)

dev.off()



###  image  ###
#outi=out[rownames(out)%in% name_image,]

# table(out[,"R"])
indexi=which(out[,"R"]<=10 & rownames(out)%in% name_image)

outij=out[indexi,]
nrow(outij)

indexi=indexi[order(outij[,"R"])]

outij=outij[order(outij[,"R"]),]
rownames(outij)

betai=list(); i=1
for (i in 1:length(indexi)) {
  betai[[i]]=sapply(BetaBoot,function(x){x[indexi[i],]})
}

# 95% CI
CIup=NULL; CIlow=NULL
for (i in 1:length(indexi)) {
  CIup=cbind(CIup,apply(betai[[i]],1,quantile,prob=0.975))
  CIlow=cbind(CIlow,apply(betai[[i]],1,quantile,prob=0.025))
}


#namei=c("ThalamusR","PallidumR","ThalamusL","PallidumL","PutamenR","AmygdalaL","CerebellumR","AccumbensR","CaudateL","AmygdalaR")
namei=c("ThalamusR","PallidumR","ThalamusL","PallidumL","PutamenR","AmygdalaL","AccumbensR","CerebellumCR")#,"CaudateL")
ii=sign(apply(outij[,1:nw0],1,mean))==1
indexi1=which(ii); indexi2=which(!ii)
coli=ifelse(ii,"black","red2");table(coli)
#coli[coli=="black"]=c("gray0","gray20","gray40","gray60")
#coli[coli=="red"]=c("red4","red3","orangered3","orangered","orange4")
ltyi=numeric(nrow(outij))
ltyi[indexi1]=1;ltyi[indexi2]=2


i=1

#pdf("Beta_image2.pdf")
#tiff(paste("Figure/Beta_image",i,".tif",sep=""),width=6,height=6,units='in',res=300)
#tiff(paste("Beta_image",h,".tif",sep=""),width=6,height=6,units='in',res=300)

cairo_ps(paste("Figure_revision/Beta_image",i,".eps",sep=""))

plot(w0,outij[indexi1[i],1:nw0],ylim=c(-2,2),col=coli[indexi1[i]],type="l",ylab="Beta",xlab=expression("CAP"[s]),lty=ltyi[indexi1[i]],main=paste("Imaging biomarkers",i),cex.main=1.7,cex.lab=1.3)
abline(h=0,col="grey")
points(w0,outij[indexi1[i],1:nw0],type="l",lty=ltyi[indexi1[i]],col=coli[indexi1[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi1[i]],rev(CIup[,indexi1[i]])),col="#D4D4D488",border=NA)
points(w0,outij[indexi2[i],1:nw0],type="l",lty=ltyi[indexi2[i]],col=coli[indexi2[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi2[i]],rev(CIup[,indexi2[i]])),col="#D4D4D488",border=NA)
abline(v=0.67,lty="dashed",col="gray30",lwd=1.5);abline(v=0.85,lty="dashed",col="gray30",lwd=1.5)

legend("topleft",legend=namei[c(indexi1[i],indexi2[i])],col=coli[c(indexi1[i],indexi2[i])],text.col=coli[c(indexi1[i],indexi2[i])],lty=c(1,2),lwd=2,cex=1.3)


dev.off()



i=2

#pdf("Beta_image2.pdf")
#tiff(paste("Figure/Beta_image",i,".tif",sep=""),width=6,height=6,units='in',res=300)
#tiff(paste("Beta_image",h,".tif",sep=""),width=6,height=6,units='in',res=300)

cairo_ps(paste("Figure_revision/Beta_image",i,".eps",sep=""))

plot(w0,outij[indexi1[i],1:nw0],ylim=c(-2,2),col=coli[indexi1[i]],type="l",ylab="Beta",xlab=expression("CAP"[s]),lty=ltyi[indexi1[i]],main=paste("Imaging biomarkers",i),cex.main=1.7,cex.lab=1.3)
abline(h=0,col="grey")
points(w0,outij[indexi1[i],1:nw0],type="l",lty=ltyi[indexi1[i]],col=coli[indexi1[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi1[i]],rev(CIup[,indexi1[i]])),col="#D4D4D488",border=NA)
points(w0,outij[indexi2[i],1:nw0],type="l",lty=ltyi[indexi2[i]],col=coli[indexi2[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi2[i]],rev(CIup[,indexi2[i]])),col="#D4D4D488",border=NA)
abline(v=0.67,lty="dashed",col="gray30",lwd=1.5);abline(v=0.85,lty="dashed",col="gray30",lwd=1.5)

legend("topleft",legend=namei[c(indexi1[i],indexi2[i])],col=coli[c(indexi1[i],indexi2[i])],text.col=coli[c(indexi1[i],indexi2[i])],lty=c(1,2),lwd=2,cex=1.3)

dev.off()


i=3

#pdf("Beta_image2.pdf")
#tiff(paste("Figure/Beta_image",i,".tif",sep=""),width=6,height=6,units='in',res=300)
#tiff(paste("Beta_image",h,".tif",sep=""),width=6,height=6,units='in',res=300)

cairo_ps(paste("Figure_revision/Beta_image",i,".eps",sep=""))

plot(w0,outij[indexi1[i],1:nw0],ylim=c(-2,2),col=coli[indexi1[i]],type="l",ylab="Beta",xlab=expression("CAP"[s]),lty=ltyi[indexi1[i]],main=paste("Imaging biomarkers",i),cex.main=1.7,cex.lab=1.3)
abline(h=0,col="grey")
points(w0,outij[indexi1[i],1:nw0],type="l",lty=ltyi[indexi1[i]],col=coli[indexi1[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi1[i]],rev(CIup[,indexi1[i]])),col="#D4D4D488",border=NA)
points(w0,outij[indexi2[i],1:nw0],type="l",lty=ltyi[indexi2[i]],col=coli[indexi2[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi2[i]],rev(CIup[,indexi2[i]])),col="#D4D4D488",border=NA)
abline(v=0.67,lty="dashed",col="gray30",lwd=1.5);abline(v=0.85,lty="dashed",col="gray30",lwd=1.5)

legend("topleft",legend=namei[c(indexi1[i],indexi2[i])],col=coli[c(indexi1[i],indexi2[i])],text.col=coli[c(indexi1[i],indexi2[i])],lty=c(1,2),lwd=2,cex=1.3)

dev.off()


i=4

#pdf("Beta_image2.pdf")
#tiff(paste("Figure/Beta_image",i,".tif",sep=""),width=6,height=6,units='in',res=300)
#tiff(paste("Beta_image",h,".tif",sep=""),width=6,height=6,units='in',res=300)

cairo_ps(paste("Figure_revision/Beta_image",i,".eps",sep=""))

plot(w0,outij[indexi1[i],1:nw0],ylim=c(-2,2),col=coli[indexi1[i]],type="l",ylab="Beta",xlab=expression("CAP"[s]),lty=ltyi[indexi1[i]],main=paste("Imaging biomarkers",i),cex.main=1.7,cex.lab=1.3)
abline(h=0,col="grey")
points(w0,outij[indexi1[i],1:nw0],type="l",lty=ltyi[indexi1[i]],col=coli[indexi1[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi1[i]],rev(CIup[,indexi1[i]])),col="#D4D4D488",border=NA)
points(w0,outij[indexi2[i],1:nw0],type="l",lty=ltyi[indexi2[i]],col=coli[indexi2[i]],lwd=2.5)
polygon(c(w0,rev(w0)),c(CIlow[,indexi2[i]],rev(CIup[,indexi2[i]])),col="#D4D4D488",border=NA)
abline(v=0.67,lty="dashed",col="gray30",lwd=1.5);abline(v=0.85,lty="dashed",col="gray30",lwd=1.5)

legend("topleft",legend=namei[c(indexi1[i],indexi2[i])],col=coli[c(indexi1[i],indexi2[i])],text.col=coli[c(indexi1[i],indexi2[i])],lty=c(1,2),lwd=2,cex=1.3)

dev.off()




###  Plot network  ###
library(igraph)

wa=c(0.67,0.85,1.0)
M0=which(w0%in%wa)
M0=c(18,36,51);w0[M0]

### image name
# namei=name_image
namei=c("Lateral","Inf","White","Cortex","Thalamus","Caudate","Putamen","Pallidum","Hippocampus","Amygdala","Accumbens","VentralDC","Vessel","Choroid")
namei=paste(rep(namei,2),rep(c("L","R"),each=length(namei)),sep="")

r=0.7;nn=1 # for ii=1
temm=apply(abs(Betai[9:p,]),1,mean)!=0
#name_image[temm]

ii=1
betaj=Betai[9:p,M0[ii]]
indexj=which(betaj!=0 & temm)
betaj=betaj[indexj]

Adjm=Omega[indexj+8,indexj+8]
Adjm=abs(Adjm);diag(Adjm)=0
Adjm[Adjm<=r]=0;
tem=Adjm;tem[tem>r]=1
#Adjm[Adjm>=r]=1
nadj=apply(tem,1,sum)

colnames(Adjm)=rownames(Adjm)=namei[indexj]

indexj=(nadj>=nn)
Adjm=Adjm[indexj,indexj]

net=graph.adjacency(Adjm,mode="undirected",weighted=T,diag=FALSE) 

V(net)$size=abs(betaj[indexj])*20
#V(net)$size=rep(30,length(betaj[indexj]))
V(net)$color=ifelse(betaj[indexj]>0,"blue","red")
E(net)$weight=rep(0.8,length(E(net)$weight))*3

# set.seed(1213)
# #pdf(paste("Network2_w",ii,".pdf",sep=""))
# cairo_ps(paste("Figure/Network_w",ii,".eps",sep=""))
# #tiff(paste("Network2_w",ii,".tif",sep=""),width=6,height=6,units='in',res=300)
# plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(net)$weight,main=list(expression(paste("CAP"[s],"=0.67",sep="")),cex=1.5),vertex.label.dist=1,vertex.label.cex=1.5)#,vertex.label.degree=pi/1.2)
# dev.off()

set.seed(1213)
#pdf(paste("Network2_w",ii,".pdf",sep=""))
cairo_ps(paste("Figure/Network_w",ii,".eps",sep=""))
#tiff(paste("Network2_w",ii,".tif",sep=""),width=6,height=6,units='in',res=300)
plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(net)$weight,main=list(expression(paste("CAP"[s],"=0.67",sep="")),cex=2),vertex.label.dist=1,vertex.label.cex=2,vertex.label.degree=-pi/1.2)
dev.off()

###
r=0.7;ic=0;nn=1 # for ii=2,3

ii=2;ii=3
for(ii in 2:3){
  
  betaj=Betai[9:p,M0[ii]]
  indexj=which(betaj!=0 & temm)
  betaj=betaj[indexj]
  
  Adjm=Omega[indexj+8,indexj+8]
  Adjm=abs(Adjm);diag(Adjm)=0
  Adjm[Adjm<=r]=0;
  tem=Adjm;tem[tem>r]=1
  #Adjm[Adjm>=r]=1
  nadj=apply(tem,1,sum)
  
  colnames(Adjm)=rownames(Adjm)=namei[indexj]
  
  indexj=(nadj>=nn)
  Adjm=Adjm[indexj,indexj]
  
  net=graph.adjacency(Adjm,mode="undirected",weighted=T,diag=FALSE) 
  
  V(net)$size=abs(betaj[indexj])*20
  #V(net)$size=rep(30,length(betaj[indexj]))
  V(net)$color=ifelse(betaj[indexj]>0,"blue","red")
  E(net)$weight=E(net)$weight*3
  
  set.seed(9)
  tiff(paste("Figure/Network_w",ii,".tif",sep=""),width=6,height=6,units='in',res=300)
  #pdf(paste("Network_w",ii,".pdf",sep=""))
  plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(net)$weight)
  dev.off()
}

ii=2
set.seed(1213)
#pdf(paste("Network2_w",2,".pdf",sep=""))
cairo_ps(paste("Figure/Network_w",2,".eps",sep=""))
plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(net)$weight,main=list(expression(paste("CAP"[s],"=0.85",sep="")),cex=2),vertex.label.dist=1,vertex.label.cex=2,vertex.label.degree=-pi/1.2)
dev.off()

ii=3
set.seed(1213)
#pdf(paste("Network2_w",3,".pdf",sep=""))
cairo_ps(paste("Figure/Network_w",3,".eps",sep=""))
plot.igraph(net,vertex.label=V(net)$name,layout=layout.fruchterman.reingold, edge.color="black",edge.width=E(net)$weight,main=list(expression(paste("CAP"[s],"=1.00",sep="")),cex=2),vertex.label.dist=1,vertex.label.cex=2,vertex.label.degree=-pi/1.2)
dev.off()





######################
#####  Analysis  #####
######################


Beta=list();Betat=list();cuts=list()


###############
###  Local  ###
###############
# range(w);quantile(w)
# c(0.67,0.85,1.0)
# w0=c(seq(0.7,0.98,by=0.02),c(0.67,1.0));w0=unique(w0);w0=w0[order(w0)]
w0=seq(0.50,1.00,by=0.01);w0=unique(w0);w0=w0[order(w0)]
hs=c(0.1,0.15,0.2);nh=length(hs)

###################
###  AdalocNet  ###
###################
hs=c(0.125,0.15,0.2);nh=length(hs)
indexh=ih=1;ia=1

aini=locoxini(x=xstd,y,w,w0,hs[ih],alambda,nalambda,rlambda=NULL,alinear=FALSE)

modeli=1;alphas=seq(0.1,0.9,by=0.1);lambda=NULL
outj=ttCoxnet(x=xstd,y=y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid,keep.beta=T)

### Trimming
cv.fiti0=outj
indexi=which.max(cv.fiti0$fit$cvm)
indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))

betati=list()
ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)

#Trim1
cutoffj=seq(0,max(absmi)-0.001,length=20)
cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid)
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
cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,adaptive=T,aini=aini,nfolds=nfolds,foldid=foldid)
# plot(cuti$cutoff,cuti$cvm,type="b")
ddi0=15
repeat{
  if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
  if(cuti$cvm[ddi0]>=cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
  ddi0=ddi0-1
}
cutoffi1=cuti$cutoff[ddi0]
ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0

Betat[[modeli]]=betati
cuts[[modeli]]=data.frame(cut1=cutoffi1)


### rank
rtem=apply(abs(Betat[[1]][[1]]),1,mean)
rtem=rank(exp(-rtem))



################
###  locNet  ###
################
hs=c(0.1,0.15,0.2);nh=length(hs)

modeli=2;alphas=seq(0.1,0.9,by=0.1);lambda=NULL

indexh=ih=2;ia=1
outj=ttCoxnet(x=xstd,y=y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)

### Trimming
cv.fiti0=outj
indexi=which.max(cv.fiti0$fit$cvm)
indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))

betati=list()
ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)

#Trim1
cutoffj=seq(0,max(absmi)-0.001,length=20)
cutj=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
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
cuti=cv.cutoff.locnet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],Omega=Omega,alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
# plot(cuti$cutoff,cuti$cvm,type="b")
ddi0=15
repeat{
  if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
  if(cuti$cvm[ddi0]>cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
  ddi0=ddi0-1
}
cutoffi1=cuti$cutoff[ddi0]
ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0

Betat[[modeli]]=betati
cuts[[modeli]]=data.frame(cut1=cutoffi1)



#################
###  locEnet  ###
#################

modeli=3;alphas=seq(0.1,0.9,by=0.1);lambda=NULL

ih=indexh=2;ia=1
outj=ttCoxnet(x=xstd,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Enet",alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)

### Trimming
cv.fiti0=outj
indexi=which.max(cv.fiti0$fit$cvm)
indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))

betati=list()
ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)

#Trim1
cutoffj=seq(0,max(absmi)-0.001,length=20)
cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
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
cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],alpha=alphas[ia],lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
# plot(cuti$cutoff,cuti$cvm,type="b")
ddi0=15
repeat{
  if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
  if(cuti$cvm[ddi0]>cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
  ddi0=ddi0-1
}
cutoffi1=cuti$cutoff[ddi0]
ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0

Betat[[modeli]]=betati
cuts[[modeli]]=data.frame(cut1=cutoffi1)



###############
###  locL1  ###
###############

modeli=4;alpha=1.0;lambda=NULL

ih=indexh=2
outj=ttCoxnet(x=xstd,y,w=w,w0=w0,h=hs[ih],hnext=c(hs[ih],hs[ih]+0.1),penalty="Enet",alpha=alpha,lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)

### Trimming
cv.fiti0=outj
indexi=which.max(cv.fiti0$fit$cvm)
indexj=max(which(cv.fiti0$fit$cvm>cv.fiti0$fit$cvm[indexi]-cv.fiti0$fit$cvse[indexi]))

betati=list()
ib=1;betati[[ib]]=cv.fiti0$Beta[[indexi]]
absmi=apply(abs(cv.fiti0$Beta[[indexi]]),1,mean)

#Trim1
cutoffj=seq(0,max(absmi)-0.001,length=20)  
cutj=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffj,nfolds=nfolds,foldid=foldid)
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
cuti=cv.cutoff.locenet(cv.fiti0$Beta[[indexi]],x=xstd,y,w,w0,h=hs[indexh],alpha=alpha,lambda=cv.fiti0$fit$lambda[indexi],cutoff=cutoffi,nfolds=nfolds,foldid=foldid)
# plot(cuti$cutoff,cuti$cvm,type="b")
ddi0=15
repeat{
  if(ddi0<=2){ddi0=which.max(cuti$cvm);break}
  if(cuti$cvm[ddi0]>cuti$cvm[ddi0+1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-1] & cuti$cvm[ddi0]>=cuti$cvm[ddi0-2])break
  ddi0=ddi0-1
}
cutoffi1=cuti$cutoff[ddi0]
ib=2;betati[[ib]]=cv.fiti0$Beta[[indexi]];betati[[ib]][absmi<=cutoffi1,]=0

Betat[[modeli]]=betati
cuts[[modeli]]=data.frame(cut1=cutoffi1)


temi=Betat[[4]][[1]]
rownames(temi)=c(namex,name_image)

########################
###  Output - local  ###
########################

wa=c(0.67,0.85,1.0)
M0=which(w0%in%wa)
M0=c(18,36,51);w0[M0]

ic=2

###  Estimates  ###
out=NULL
for(k in 1:4){out=cbind(out,NA,Betat[[k]][[ic]][,M0])}
out[out==0]=NA

rtem=apply(abs(Betat[[1]][[ic]]),1,mean)
rtem=rank(exp(-rtem))
out=cbind(out,NA,rtem)

rownames(out)=c(namex,name_image)
write.csv(out,"R_tmp.csv",na="")


###  non-zero  ###
nzero=NULL
for(k in 1:4){
  nzeroi=apply(Betat[[k]][[ic]][,M0],2,function(x){sum(x!=0)})
  nzero=c(nzero,NA,nzeroi)
}

write.csv(t(nzero),"R_tmp.csv",na="")


###  C index  ###
# https://www.lerner.ccf.org/qhs/outcomes/documents/Gonen.pdf
nn=50;hi=0.02

xstdi=xstd;yi=y
betai=Betat[[1]][[ic]][,1]
scorei=xstdi%*%betai

ib=iw=1
result=NULL;inc=NULL
for(ib in 1:length(Betat)){
  resulti=NULL
  for(iw in M0){
    
    indexi=abs(w-w0[iw])<=hi
    xstdi=xstd[indexi,];yi=y[indexi,]
    betai=Betat[[ib]][[ic]][,iw]
    scorei=xstdi%*%betai
    
    den=0;num=0
    for(i in 1:(nrow(xstdi)-1)){
      for(j in (i+1):nrow(xstdi)){
        Iij=(yi$time[i]<yi$time[j] & yi$status[i]==1)
        Iji=(yi$time[j]<yi$time[i] & yi$status[j]==1)
        if(Iij+Iji==1){
          den=den+1
          num=num+Iij*(scorei[i]>scorei[j])+Iji*(scorei[j]>scorei[i])
        }
      }
    }
    inc=c(inc,den)
    resulti=c(resulti,num/den)
  }
  result=c(result,NA,resulti)
}
result
inc

write.csv(t(result),"R_tmp.csv",na="")



###################
###  non-Local  ###
###################

#####  initial for AdaNet  #####
outi=ttCoxnet(xstd,y,penalty="Enet",alpha=0.0,lambda=0.0) # MLE
aini2=list(wbeta=1/abs(outi$Beta),sgn=sign(outi$Beta))

# library(survival)
# outi=ttCoxnet(xstd,y,penalty="Enet",alpha=0.0,lambda=0.0)
# cfit2=coxph(Surv(y$time,y$status)~xstd)
# cbind(outi$Beta,cfit2$coef)
# all.equal(outi$Beta,cfit2$coef);range(outi$Beta-cfit2$coef)
# name_image[which(summary(cfit2)$coefficients[,"Pr(>|z|)"]<0.05)]
# write.csv(summary(cfit2)$coefficients,"R_tmp.csv")
# # "Left_Cerebellum_Cortex";"Left_Pallidum";"Left_Amygdala";"Right_Inf_Lat_Vent";"Right_Cerebellum_White_Matter"


#####  Estimation  #####
Beta2=list();select.para2=list()

###  ANet  ###
modeli=1;alphas=seq(0.1,0.9,by=0.1);lambda=NULL

# cll=makeCluster(ncluster);registerDoSNOW(cll)
outi=list()
for(ia in 1:nalpha) {
  fiti=ttCoxnet(xstd,y,penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,adaptive=T,aini=aini2,nfolds=nfolds,foldid=foldid,keep.beta=T)
  outi[[ia]]=(fiti)
}
# stopCluster(cll)

### Max
indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
indexi=which.max(outi[[indexa]]$fit$cvm)

temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli)
select.para2[[modeli]]=temparai

Beta2[[modeli]]=outi[[indexa]]$Beta[,indexi]


###  Net  ###
modeli=2;alphas=seq(0.1,0.9,by=0.1);lambda=NULL

cll=makeCluster(ncluster);registerDoSNOW(cll)
outi=list()
for(ia in 1:nalpha) {
  fiti=ttCoxnet(xstd,y,penalty="Network",Omega=Omega,alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
  outi[[ia]]=(fiti)
}
stopCluster(cll)

### Max
indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
indexi=which.max(outi[[indexa]]$fit$cvm)

temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli)
select.para2[[modeli]]=temparai

Beta2[[modeli]]=outi[[indexa]]$Beta[,indexi]


###  Enet  ###
#glmi=glmnet(x,y,family="cox",alpha=alpha,lambda=lambda)
modeli=3;alphas=seq(0.1,0.9,by=0.1);lambda=NULL

# cll=makeCluster(ncluster);registerDoSNOW(cll)
outi=list()
for(ia in 1:nalpha) {
  fiti=ttCoxnet(xstd,y,penalty="Enet",alpha=alphas[ia],lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)
  outi[[ia]]=(fiti)
}
# stopCluster(cll)

### Max
indexa=which.max(sapply(outi,function(x){max(x$fit$cvm)}))
indexi=which.max(outi[[indexa]]$fit$cvm)
indexij=which(outi[[indexa]]$fit$index=="se")

temparai=data.frame(alpha=alphas[indexa],lambda=outi[[indexa]]$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli)
select.para2[[modeli]]=temparai

Beta2[[modeli]]=outi[[indexa]]$Beta[,indexi]


###  Lasso  ###
modeli=4;alpha=1.0;lambda=NULL

fiti=ttCoxnet(xstd,y,penalty="Enet",alpha=alpha,lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,keep.beta=T)

indexi=which.max(fiti$fit$cvm)

temparai=data.frame(alpha=alpha,lambda=fiti$fit$lambda[indexi],h=NA,cvh0=NA,cvh1=NA,model=modeli)
select.para2[[modeli]]=temparai

Beta2[[modeli]]=outi[[indexa]]$Beta[,indexi]


###  None  ###
# library(survival)
modeli=5;alpha=1.0;lambda=0.0

fiti=ttCoxnet(xstd,y,penalty="Enet",alpha=alpha,lambda=lambda)
Beta2[[modeli]]=fiti$Beta


###  Output  ###
out=do.call("cbind",Beta2)
rownames(out)=colnames(x)
colnames(out)=c("ANet","Net","Enet","Lasso","None")

nz=sapply(Beta2,function(x){sum(x!=0)})

## C-index
result=NULL;inc=NULL
for(i in 1:5){
  xstdi=xstd;yi=y
  betai=Beta2[[i]]
  scorei=xstdi%*%betai
  
  den=0;num=0
  for(i in 1:(nrow(xstdi)-1)){
    for(j in (i+1):nrow(xstdi)){
      Iij=(yi$time[i]<yi$time[j] & yi$status[i]==1)
      Iji=(yi$time[j]<yi$time[i] & yi$status[j]==1)
      if(Iij+Iji==1){
        den=den+1
        num=num+Iij*(scorei[i]>scorei[j])+Iji*(scorei[j]>scorei[i])
      }
    }
  }
  inc=c(inc,den)
  result=c(result,num/den)
}

out[out==0]=NA

write.csv(rbind(out,nz,result),"R_tmp.csv",na="")




