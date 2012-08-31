#############################################
### Functions for testing significance of ###
### per-gene categorized mutation rates   ###
#############################################

# Fetch command line arguments
input_file = as.character(commandArgs()[4]);
output_file = as.character(commandArgs()[5]);
run_type = as.character(commandArgs()[6]);

gethist=function(xmax,n,p,ptype="positive_log")
{
  dbinom(0:xmax,n,p)->ps
  ps=ps[ps>0]
  lastp=1-sum(ps)
  if (lastp>0) ps=c(ps,lastp)
  if (ptype=="positive_log") ps=-log(ps)
  ps
}

binit=function(x,hmax,bin,dropbin=T)
{
  bs=as.integer(x/bin)
  bs[bs>hmax/bin]=hmax/bin
  bs[is.na(bs)]=hmax/bin
  tapply(exp(-x),as.factor(bs),sum)->bs
  bs=bs[bs>0]
  bs=-log(bs)
  if (dropbin) bs=as.numeric(bs)
  bs
}

convolute_b=function(a,b)
{
  tt=NULL
  for (i in a){
    for (j in b){
      temp=i+j
      tt=c(tt,temp)
  }}
  tt
}

mut_class_test=function(x,xmax=100,hmax=25,bin=0.001)
{
  x=as.data.frame(x)
  colnames(x)=c("n","x","e")
  x$p=NA
  x$lh0=NA
  x$lh1=NA
  hists=NULL
  for (i in 1:nrow(x))
  {
    x$p[i]=binom.test(x$x[i],x$n[i],x$e[i],alternative="greater")$p.value
    x$lh0[i]=dbinom(x$x[i],x$n[i],x$e[i],log=T)
    x$lh1[i]=dbinom(x$x[i],x$n[i],x$x[i]/x$n[i],log=T)
    ni=x$n[i];ei=x$e[i]
    gethist(xmax,ni,ei,ptype="positive_log")->bi
    binit(bi,hmax,bin)->bi
    if (i==1) hist0=bi
    if (i>1 & i<nrow(x)) {hist0=convolute_b(hist0,bi);binit(hist0,hmax,bin)->hist0}
    if (i==nrow(x)) hist0=convolute_b(hist0,bi)
  }

  # Fisher combined p-value
  q= (-2)*sum(log(x$p))
  df=2*length(x$p)
  p.fisher= 1-pchisq(q, df)

  # Likelihood ratio test
  q=2*(sum(x$lh1)-sum(x$lh0))
  df=sum(x$lh1!=0)
  if (df>0) p.lr= 1-pchisq(q, df)
  if (df==0) p.lr=1

  # Convolution test
  tx=sum(x[,"x"])
  tn=sum(x[,"n"])
  (bx=-sum(x[,"lh0"]))
  (p.convol=sum(exp(-hist0[hist0>=bx])))
  (qc=sum(exp(-hist0)))

  # Return results
  rst=list(hists=hist0,x=cbind(x,tn,tx,p.fisher,p.lr,p.convol,qc))
  rst
}

smg_test=function(in.file,test.file)
{
  source("/gscuser/qzhang/gstat/mut_class_test.R");

  read.table(in.file,header=T,sep="\t")->mut
  mut$BMR=as.numeric(as.character(mut$BMR))

  #select the rows with BMR data
  mut=mut[mut$BMR>0 & !is.na(mut$BMR) & mut$Bases>0,]
  tt=NULL
  for (gene in unique(as.character(mut$Gene)))
  {
    mutgi=mut[mut$Gene==gene,]
    mut_class_test(mutgi[,3:5],hmax=25,bin=0.001)->z
    tt=rbind(tt,cbind(mutgi,z$x[,-(1:3)]))
  }
  write.table(tt,file=test.file,quote=FALSE,row.names=F,sep="\t")
}

smg_fdr=function(in.file,fdr.file)
{
  read.table(in.file,header=T,sep="\t")->x
  x=unique(x)

  #Calculate FDR measure and write FDR output
  p.adjust(x[,2],method="BH")->fdr.fisher
  p.adjust(x[,3],method="BH")->fdr.lr
  p.adjust(x[,4],method="BH")->fdr.convol
  x=cbind(x,fdr.fisher,fdr.lr,fdr.convol)
  write.table(x,file=fdr.file,sep="\t",quote=FALSE,row.names=F)
}
