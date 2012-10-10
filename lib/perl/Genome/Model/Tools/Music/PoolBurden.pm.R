################################################
### Functions for running Pooled Burden Test ###
###############################################

# Fetch command line arguments
args = commandArgs(trailingOnly = TRUE);
input_file = as.character(args[1]);
output_file = as.character(args[2]);
processors = as.numeric(args[3]);


# See if we have the necessary packages installed to run in parallel
is.installed <- function( mypkg ) is.element( mypkg, installed.packages()[,1] );
parallel = FALSE;
if( processors > 1 & is.installed( 'doMC' ) & is.installed( 'foreach' )) {
    parallel = TRUE;
}



# hypergeometric distribution
p.hyper=function(a,b,c,d,log=F) {
  n = a+b+c+d;
  pi = lgamma(a+b+1)+lgamma(c+d+1)+lgamma(a+c+1)+lgamma(b+d+1);
  pi = pi-lgamma(a+1)-lgamma(b+1)-lgamma(c+1)-lgamma(d+1)-lgamma(n+1);
  if (!log) pi=exp(pi);
  pi; 
}

###
convolute_b=function(a,b) {
  tt=NULL;
  for (j in b) tt=c(tt,(a+j));
  tt;
}

pool.RV.test=function(x,exact.lim=1000000,simu.n4p=10000000) {

  #x=zi
  x=as.data.frame(x)
  n4.col=c("x0","x1","n0","n1")
  x$p=NA
  x$lh0=NA
  x$lh1=NA
  x$dhyper=NA
  x$dbinom=NA
  x$T=NA
  x$c=NA
  x$vs=NA
  hist0.exact=NA
  hist0=NULL
  hist0.exactb=NA
  hist0b=NULL
  pn=x$x0+x$x1+1;pn=exp(sum(log(pn)))

  for (i in 1:nrow(x)) {

    #for LR
    e0=(x$x0[i]+x$x1[i])/(x$n0[i]+x$n1[i])
    x$lh0[i]=dbinom(x$x0[i],x$n0[i],e0,log=T) + dbinom(x$x1[i],x$n1[i],e0,log=T)
    e0=x$x0[i]/x$n0[i]
    e1=x$x1[i]/x$n1[i]
    lh.0=dbinom(x$x0[i],x$n0[i],e0,log=T)
    lh.1=dbinom(x$x1[i],x$n1[i],e1,log=T)
    x$lh1[i]=lh.0+lh.1
    
    xi=as.numeric(x[i,n4.col])

    #fisher single variant test
    a=xi[1];b=xi[2];c=xi[3]-xi[1];d=xi[4]-xi[2]
    x$p[i]=fisher.test(matrix(c(a,b,c,d),2,2))$p.value

    #Pi for Hypergeometric Exact Prob test
    x$dhyper[i]=p.hyper(a,b,c,d,log=T)
    gethist2(xi,ptype="positive_log")->bi

    #Pi for Binomial Exact Prob test
    x$dbinom[i]=dbinom(a,a+b,(a+c)/(a+b+c+d),log=T)
    -dbinom(0:(a+b),a+b,(a+c)/(a+b+c+d),log=T)->bib

    #exact distribution
    if (pn < exact.lim) {
      if (i==1) hist0.exact=bi
      if (i>1) hist0.exact=convolute_b(hist0.exact,bi)
      if (i==1) hist0.exactb=bib
      if (i>1) hist0.exactb=convolute_b(hist0.exactb,bib)
    }

    #for simulated empirical distribution
    hist0=c(hist0,list(bi))
    hist0b=c(hist0b,list(bib))

    #for C-alpha test
    yi=xi[1];ni=xi[1]+xi[2];p0=xi[3]/(xi[3]+xi[4])
    x$T[i]=(yi-ni*p0)^2-ni*p0*(1-p0)
    u=c(0:ni)
    x$c[i]= sum(((u-ni*p0)^2-ni*p0*(1-p0))^2*dbinom(u,ni,p0))

    #for CAST
    (xi[1]/xi[3])*(1-xi[1]/xi[3])*xi[3] -> vs1
    (xi[2]/xi[4])*(1-xi[2]/xi[4])*xi[4] -> vs2
    xi[3]/xi[4]-> c12
    x$vs[i]=vs1+vs2*c12*c12

    #binit(bi,bin,hmax)->bi
    #if (i==1) hist0=bi
    #if (i>1 & i<nrow(x)) {hist0=convolute_b(hist0,bi);binit(hist0,bin,hmax)->hist0}
    #if (i==nrow(x)) hist0=convolute_b(hist0,bi)
  }

  # bonferroni
  p.bon=as.numeric(min(x$p)*nrow(x)); if (p.bon>1) p.bon=1

  #sum test
  a=sum(x$x0);b=sum(x$x1)
  sumi=c(a,b,x$n0[1]-a,x$n1[1]-b)
  (p.sum=fisher.test(matrix(as.numeric(sumi),2,2))$p.value)

  #exclusive sum test
  a=sum(x$x0[x$x1==0]);b=sum(x$x1[x$x0==0])
  sumi=c(a,b,x$n0[1]-a,x$n1[1]-b)
  (p.esum=fisher.test(matrix(as.numeric(sumi),2,2))$p.value)

  #CAST
  S=sum(x$x0-x$x1*x$n0/x$n1)
  SE=sqrt(sum(x$vs))
  q=S/SE
  (p.CAST.greater= pnorm(q))
  q=abs(S/SE)
  (p.CAST= 2*pnorm(-q))


  #C-alpha test
  q=sum(x$T)/sqrt(sum(x$c))
  (p.Calpha= pnorm(-q))

  #fisher combined p-value
  q=(-2)*sum(log(x$p))
  df=2*length(x$p)
  (p.fisher= 1-pchisq(q, df/2))

  #likelihood ratio test
  q=2*(sum(x$lh1)-sum(x$lh0))
  (p.lr= 1-pchisq(q, df/2))

  #logit combined p-value test
  p0=x$p;p0=p0[!is.na(p0) & p0>0 & p0<1]
  q=mean(log(p0/(1-p0)))
  df=length(p0)
  (p.logit= pt(q, df))

  #sum chisq combined p-value test
  p0=x$p;p0=p0[!is.na(p0)]
  q=sum(qchisq(p0,1))
  df=length(p0)
  (p.sumchi= pchisq(q, df))

  ###Exact prob. test H-EPT
  (bx=-sum(x$dhyper));
  p.exact=NA; qc.exact=NA;
  if (!is.na(hist0.exact)[1]) {
    p.exact=1-sum(exp(-hist0.exact[hist0.exact<bx]));
    if (p.exact<0) p.exact=sum(exp(-hist0.exact[hist0.exact>bx]));
    qc.exact=sum(exp(-hist0.exact))
  }

  ###Exact prob. test B-EPT
  (bxb=-sum(x$dbinom))
  p.exactb=NA; qc.exactb=NA
  if (!is.na(hist0.exactb)[1]) {
    p.exactb=1-sum(exp(-hist0.exactb[hist0.exactb<bxb]))
    if (p.exactb<0) p.exactb=sum(exp(-hist0.exactb[hist0.exactb>bxb]))
    qc.exactb=sum(exp(-hist0.exactb))
  }

  ###simulation test sH-EPT
  ds=0;di=0
  for (i in hist0) {
    sample(i,simu.n4p,replace=T,prob=exp(-i))->di
    ds=ds+di
  }
  p.empi=sum(ds>=bx)/length(ds)


  ###simulation test sB-EPT
  ds=0;di=0
  for (i in hist0b) {
    sample(i,simu.n4p,replace=T,prob=exp(-i))->di
    ds=ds+di
  }
  p.empib=sum(ds>=bxb)/length(ds)

  X0=sum(x$x0)
  X1=sum(x$x1)
  N0=mean(x$n0)
  N1=mean(x$n1)
  MAF0=X0/N0
  MAF1=X1/N1
  VN=nrow(x)
  
  # return
  rst=list(x=x,
    p=cbind(VN,X0,X1,N0,N1,MAF0,MAF1,pn,p.bon,p.CAST,p.CAST.greater,p.Calpha,p.sum,p.esum,p.fisher,p.logit,p.sumchi,p.lr,p.exact,p.empi,p.exactb,p.empib),
    qc=qc.exact,qcb=qc.exactb);
  rst

}


# hypergeometric (case + control)
# two binom
gethist2=function(xi,ptype="positive_log") {
  xi=as.numeric(xi);
  x0=xi[1];x1=xi[2];n0=xi[3];n1=xi[4];
  nx=x0+x1;
  ps=NULL;
  for (a in 0:nx) { 
    b=nx-a;c=n0-a;d=n1-b;
    ps=c(ps,p.hyper(a,b,c,d,log=T))
  }

  #ps=ps[ps>0]
  #lastp=1-sum(ps)
  #if (lastp>0) ps=c(ps,lastp)
  if (ptype=="positive_log") ps=-ps
  ps;
}

dotest <- function( idx, input.data, zgenes ) {
    step = round( length( zgenes ) / processors );
    start = step * ( idx - 1 ) + 1;
    stop = step * idx;
    if( idx == processors ) { stop = length( zgenes ); }
    tt = NULL;
    for( Gene in zgenes[start:stop] ) {
        mutgi = input.data[input.data$gene==Gene,];
        pi=NA;
        pi=try(pool.RV.test(mutgi,1000000,1000000)$p);
        if (class(pi)!="try.error") tt=rbind(tt,pi);
        #mut_class_test( mutgi[,2:5], hmax = 25, bin = 0.001 ) -> z;
        #tt = rbind( tt, cbind( Gene, unique( z$x[,(9:11)] )));
    }
    return( tt );
}

pool_burden_test <- function (in.file,pval_file) {

  z <- read.csv(in.file);
  z <- z[z$x0+z$x1>0,]

  #uniq_genes = a list of gene names
  tg <- table(z$gene);
  gs <- names(tg)[tg>=1];
  gs <- gs[nchar(gs)>0];
  uniq_genes <- gs[!is.na(gs)];
  
  tt = NULL;
  # Run in parallel if we have the needed packages, or fall back to the old way
  if( parallel ) {
    library( 'doMC' );
    library( 'foreach' );
    registerDoMC();
    cat( "Parallel backend installed - splitting across", processors, "cores\n" );

    options( cores = processors );
    mcoptions <- list( preschedule = TRUE );

    #zgenes = unique( as.character( mut$Gene ));
    tt = foreach( idx = 1:processors, .combine="combineresults", .options.multicore = mcoptions ) %dopar% {
      dotest( idx, z, uniq_genes);
    }

   

    #write.table( tt, file = pval_file, quote = FALSE, row.names = F, sep = "\t" );
  }
  else {
    for( Gene in uniq_genes) {
      mutgi = z[z$gene==Gene,];
      pi=NA;
      pi=try(pool.RV.test(mutgi,1000000,1000000)$p);
      if (class(pi)!="try.error") { tt=rbind(tt,pi); }
    }
    #write.table(tt, file=pval_file, na="", quote=T, sep="\t",row.names = TRUE,)
    #write.table( tt, file = pval_file, quote = FALSE, row.names = F, sep = "\t" );
  }

  rownames(tt)=gs;
  CAST=tt[,"p.CAST"];
  CAST.greater=tt[,"p.CAST.greater"];
  Calpha=tt[,"p.Calpha"];
  TFT=tt[,"p.sum"];
  ETFT=tt[,"p.esum"];
  LRT=tt[,"p.lr"];
  EPT=tt[,"p.exact"];
  sEPT=tt[,"p.empi"];

  tt=cbind(tt[,1:7],CAST,CAST.greater,Calpha,TFT,ETFT,LRT,EPT,sEPT);
  tt=signif(tt,6);

  write.csv(tt, file=pval_file, na="", quote=F)
  #write.table(tt, file=pval_file, na=" ", quote=F, sep="\t");

}


pool_burden_test(input_file,output_file);




  

