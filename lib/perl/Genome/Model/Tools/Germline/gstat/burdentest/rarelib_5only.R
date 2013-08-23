# we want
# aSum + pCMC + SPWST + SKAT + VT
 

library(SKAT)

############# funtions ############

#############################################
collapse.analysis=function(z)
{
    #################### prepare data 
    #z=samp
    #z$x=z$x[,sd(z$x,na.rm=T)>0]
    design=z$design
    id=z$id
    grp=z$grp
    #y=z$y
    #x=z$x
    ped=z$ped
    trait=z$trait
    variants=z$variants
    covar=z$covar
    covars=z$covars
    ytype=z$ytype
    #x[x>1]=1
    data=z$data
    x=data[,variants]
    y=data[,trait]
    bls = z$bls


    ######################### single variant analysis 
    svt=NULL
    xnames=colnames(x)
    if (!is.null(bls)) {
        svt <- bls.pvalues(bls, data[,trait])
    } else {
        for (xi in xnames)
        {
            #fit=try(cor(x[,xi],y,use="pairwise.complete.obs"))
            #fit=try(cor.test(x[,xi],y,alternative="less"))
            #print(xi)
            fit=try(myglm(z=data,trait,variant=xi,covar,ytype,covar_matrix=z$covar_matrix))
            if (class(fit)!="try-error") svt=rbind(svt,fit$coef)
        }
        xnames=sort(intersect(xnames,rownames(svt)))
        svt=svt[xnames,]
    }

    r1=svt[,1];p2=svt[,4];p1=p2/2;p1[r1>0]=1-p1[r1>0]
    x=t(data[,xnames])

    N=ncol(x)
    MAF=rowSums(x)/N/2
    MAC=rowSums(x)

    svt=cbind(svt,N,MAF,MAC)

    ######################## skat ,mreg #################

    ####skat & mreg
    call.skat=function(Z,Y,ytype="C",X=NULL)
    {
        ###SKAT
        Z=as.matrix(Z)
        Y=as.matrix(Y)
        md="Y~1"
        if (!is.null(X)) {
            X=as.matrix(X);md="Y~X"
        }
        obj<-SKAT_Null_Model(formula(md), out_type=ytype)
        SKAT(Z, obj,weights.beta=c(1,1))$p.value->skat
        #SKAT(Z, obj,kernel = "linear.weighted", weights.beta=c(1,25))$p.value->skat
        #SKAT(Z, obj)$p.value->skat
        #### multi var, lm, lme
        #if (ytype=="C") y=obj$res
        #if (ytype=="D") y=obj$re1$res
        #md="y"; for (i in colnames(Z)) md=paste(md,i,sep="+")
        #formula(sub("[+]","~",md))->md
        #lm(md,data=as.data.frame(cbind(y,Z)))->fit
        #summary(fit)$fstat->sfit
        #1-pf(sfit[1],sfit[2],sfit[3])->mreg
        p=list(skat=skat,mreg=NA)
        p
    }

    ### call skat & mreg
    mv=c(NA,NA)
    if ( !("is.permu" %in% names(z)) )
    {
       if (ytype=="Q") ytp="C"
       if (ytype=="B") ytp="D"
       if (!is.null(covars))  try(call.skat(Z=as.matrix(data[,xnames]),Y=as.matrix(data[,trait]),ytype=ytp,X=as.matrix(data[,covars]))) -> fit
       if (is.null(covars))  try(call.skat(Z=as.matrix(data[,xnames]),Y=as.matrix(data[,trait]),ytype=ytp,X=NULL)) -> fit
       if (class(fit)!="try-error") mv=c(fit$skat,fit$mreg)
    } #disable skat in permutation 
    names(mv)=c("skat","mreg")

    #### Price's VT ####

    call.VT.price=function(Z,Y,X=NULL)
    {
        if (!is.null(X)) lm(Y~X)$resid->Y
        fs=maf(Z);p0=1;fcut0=NA
        for (fi in sort(unique(fs)))
        {
            cor.test(rowSums(as.matrix(Z[,fs<=fi])),Y)$p.val->pi
            if (pi<p0) {p0=pi; fcut0=fi}
        }
        p=list(vt=p0,fcut=fcut0)
        p
    }
    ### call VT
    if (!is.null(covars))  try(call.VT.price(Z=as.matrix(data[,xnames]),Y=data[,trait],X=as.matrix(data[,covars]))) -> fit
    if (is.null(covars))  try(call.VT.price(Z=as.matrix(data[,xnames]),Y=data[,trait],X=NULL)) -> fit
    vt=c(NA,NA)
    if (class(fit)!="try-error") vt=c(fit$vt,fit$fcut)
    names(vt)=c("VT","VT_MAF")

    ################## weighting, collapsing, sum score ###########################

    ############################################################### for any designs

    w=1;colSums(x*w,na.rm=T)->xi; cast=NA;   #cast=xi # CAST
    xi[xi>1]=1;cmc=xi #CMC

#w=as.numeric(r1<0);colSums(x*w,na.rm=T)->xi;xi[xi>1]=1;xdw=xi 
#w=as.numeric(r1>0);colSums(x*w,na.rm=T)->xi;xi[xi>1]=1;xup=xi

#w=p1-0.5; colSums(x*w,na.rm=T)->pwst # PWST rescaled uniform
#w=log(p1/(1-p1));colSums(x*w,na.rm=T)->pwstlog  # PWST rlogit escaled normal

    w=p2; w[p2<0.1 & r1<0]=-1; w[w>0]=1; colSums(x*w,na.rm=T)->asum  # Han & Pan P<0.1 
    w=-log(p2)*as.numeric(r1<0);colSums(x*w,na.rm=T)->xdwlog   #SPWST  -log(p)
    w=-log(p2)*as.numeric(r1>0);colSums(x*w,na.rm=T)->xuplog

#w=r1; w[w<0]=-1; w[w>0]=1; colSums(x*w,na.rm=T)->wbd  # weighted only by direction + or - 
    #vSelect(x=t(x),y=y,w=w,method="hoffmann")->stepup
#if (ytype=="Q") {w=sumweights(t(x))$wi;colSums(x*w,na.rm=T)->wss}  # WSS population
#if (ytype=="B") {w=sumweights(t(x),ctrid=(y==0))$wi;colSums(x*w,na.rm=T)->wss}  # WSS case(1)-control(0)

cast = xdw = xup = pwst = pwstlog = wbd = wss =NA;

    xc=cbind(cast,cmc,pwst,pwstlog,asum,wbd,wss,xdw,xup,xdwlog,xuplog)

    ms=list(
        SUM=list(var="cast",test="gl",alt="two.sided",id=!is.na(y)),
        C01=list(var="cmc",test="glm",alt="two.sided",id=!is.na(y)), 
        xdw=list(var="xdw",test="gl",alt="two.sided",id=!is.na(y)),
        xup=list(var="xup",test="gl",alt="two.sided",id=!is.na(y)),
        PWST=list(var="pwst",test="gl",alt="two.sided",id=!is.na(y)), 
        aSum=list(var="asum",test="glm",alt="two.sided",id=!is.na(y)), 
        #PWST0=list(var="xp",test="cor",alt="two.sided",id=!is.na(y)),
        xdwlog=list(var="xdwlog",test="glm",alt="two.sided",id=!is.na(y)),
        xuplog=list(var="xuplog",test="glm",alt="two.sided",id=!is.na(y)),
        #PWST2=list(var="xplog",test="cor",alt="two.sided",id=!is.na(y)),
        #WBD=list(var="wbd",test="glm",alt="two.sided",id=!is.na(y)),
        WSS=list(var="wss",test="gl",alt="two.sided",id=!is.na(y))
    )

    #################################################################  different designs
    xc1=NULL;ms1=NULL
    if (design=="RDM") {
        xc1=NULL;ms1=NULL
    }

    if (design %in% c("LR") ) {
        #w=sumweights(t(x[,grp=="L"]))$wi;colSums(x*w,na.rm=T)->xwR
        #w=sumweights(t(x[,grp=="R"]))$wi;colSums(x*w,na.rm=T)->xwL
        #xc1=cbind(xwR,xwL)
        #ms1=list(
        #xwR=list(var="xwR",test="cor",alt="two.sided",id=!is.na(grp)), 
        #xwL=list(var="xwL",test="cor",alt="two.sided",id=!is.na(grp)) 
        #)
        #x0c=list(var="x0",test="fisher",alt="two.sided",id=!is.na(grp)),
        #xdwc=list(var="xdw",test="fisher",alt="less",id=!is.na(grp)),
        #xupc=list(var="xup",test="fisher",alt="greater",id=!is.na(grp)),
    }
    if (design=="LCR") {
        #w=sumweights(t(x[,grp=="C"]))$wi;colSums(x*w,na.rm=T)->xw
        #xc1=cbind(xw)
        #ms1=list(
        #xwR=list(var="xw",test="cor",alt="two.sided",id=(grp!="L")),
        #xwL=list(var="xw",test="cor",alt="two.sided",id=(grp!="R"))
        #)
        #x0Rc=list(var="x0",test="fisher",alt="two.sided",id=(grp!="L")),
        #x0Lc=list(var="x0",test="fisher",alt="two.sided",id=(grp!="R")),
        #xdwc=list(var="xdw",test="fisher",alt="less",id=!is.na(grp)),
        #xupc=list(var="xup",test="fisher",alt="greater",id=!is.na(grp)),
    } 
    ms=c(ms,ms1)
    xc=cbind(xc,xc1)

    ####################################### trait and sum score association
    pp=NULL
    for (i in names(ms))
    {
        vi=ms[[i]]$var
        if (vi %in% colnames(xc)) {
            pi=NA
            ms[[i]]$id->ids
            ms[[i]]$alt->alti
            fit=NA

            if (ms[[i]]$test=="cor" & is.null(ped) ) fit=try(cor.test(xc[ids,vi],y[ids],alternative=alti))
            if (ms[[i]]$test=="cor" & !(is.null(ped)) ) fit=try(kin.test(x=xc[ids,vi],y=y[ids],id=id[ids],ped=ped))
            if (ms[[i]]$test=="fisher") fit=try(fisher.test(table(xc[ids,vi],grp[ids]),alternative=alti))
            if (ms[[i]]$test=="glm") fit=try(myglm(cbind(data,xc)[ids,],trait=trait,variant=vi,covar=covar,ytype=ytype,covar_matrix=z$covar_matrix))

if (!is.na(fit))
{
             if (class(fit)!="try-error") {
                pi=fit$p.value
                if (!is.null(ped)) ped$theta=fit$theta
            }
}

            names(pi)=i
            pp=c(pp,pi)
        }
    }

    pmin=min(pp["xdw"],pp["xup"],na.rm=T);names(pmin)="SWBD"
    pminlog=min(pp["xdwlog"],pp["xuplog"],na.rm=T);names(pminlog)="SPWST"
    #StepUp=-stepup$chi;names(StepUp)="StepUp"
    #StepUp.PW=-stepup.pw$chi;names(StepUp.PW)="StepUp.PW"

    pp=c(pp,pmin,pminlog)
    pp=c(pp,mv,vt) 

    rst=list(z=z,svt=svt,r1=r1,p1=p1,p2=p2,p=pp)

    if (!is.null(ped)) rst$theta=ped$theta

    rst
}

################################

inv_sample <- function(x) {
    order(x)
}

collapse.test=function(z,permu=0)
{
    if (!is.null(z$covar)) {
        covar_matrix <- model.frame(formula=formula(paste(z$trait,"~",z$covar,"-1")), z$data);
        mt <- attr(covar_matrix, "terms") # allow model.frame to have updated it
        z$covar_matrix <- model.matrix(mt, covar_matrix, NULL);
    }

    collapse.analysis(z)->rst

    if (permu>0) 
    {
	z$is.permu=T
        if (!is.null(z$ped)) z$ped$theta=rst$theta
        p=rst$p
        pp=p*0
        ppp=NULL

        # compute permutations in advance
        print("Precomputing permutations...")
        permutations = mat.or.vec(permu, length(z$allid))
        tm_start <- proc.time()
        for (i in 1:permu) {
            permutations[i,] <- sample(1:length(z$allid))
        }
        tm_end <- proc.time()
        print("Precomputed permutations in (seconds):")
        print(tm_end - tm_start)

        procs = min(num.cores, permu)
        perms_per_process = as.integer(permu / procs)
        tm_start <- proc.time()
        mcopts <- options(set.seed = TRUE)
        ppp <- foreach (i = 1:procs, .combine = rbind) %dopar% {
            pmatrix = NULL;
            start = (i-1)*perms_per_process + 1
            end = i*perms_per_process
            if (i == procs)
                end = permu
            print(paste("Process ", i, " (pid ", Sys.getpid(),
                ") running permutations (",start,"-",end,")", sep=""))
            
            for (j in 1:end) {
                rownames(z$data)=z$id
                id=z$allid[permutations[j,]]
                id=id[id %in% z$id]
                if (is.null(z$covar)) {
                    z$data[,z$trait]=z$data[id,z$trait]
                } else {
                    z$data[,z$variants]=z$data[id,z$variants]
                }
                if (j < start)
                    next;

                collapse.analysis(z)$p->pi
                pi[is.na(pi)]=1
                pmatrix = rbind(pmatrix, pi)
            }
            pmatrix
        }
        for (i in 1:permu) {
            pp = pp + as.integer(ppp[i,] <= p);
        }

        rst$pp=pp/permu
        rst$ppp=ppp

        tm_end <- proc.time()
        print(paste("Processed", permu, " permutations in (seconds):"))
        print(tm_end - tm_start)

    }

    rst

}

###### kinship model #######

kin.test=function(x,y,id,ped)
{
### DATA INPUT ###

#load("kin.rdata")

sub.ID=ped$subid
mom.ID=ped$momid
dad.ID=ped$dadid
ped$kmat->kmat

x=as.data.frame(cbind(id,x,y));colnames(x)[1]=sub.ID
x=x[(x[,sub.ID] %in% colnames(kmat)),]

#library("kinship")
##### generate family id ######
#cfam <- makefamid(p[,sub.ID], p[,mom.ID], p[,dad.ID])
#cbind(cfam,p)->p
##### kinship matrix #####
#kmat <- makekinship(p$cfam,p[,sub.ID], p[,mom.ID], p[,dad.ID])

pp=NULL
for (yi in "y") {
for (xi in "x") {
markname=xi
phenotype=yi
nxy=sum( (!is.na(x[,xi])) & (!is.na(x[,yi])))
Beta=NA;SE=NA;t=NA;p.value=NA;DF=NA;sampleN=NA;Pval0=NA;rsquare=NA;theta=ped$theta
if (nxy>10) 
{
fix.eff=paste(yi,"~",xi)
#if (!is.null(covariates)) {for (covi in covariates) fix.eff=paste(fix.eff,"+",covi) }
fix.eff=formula(fix.eff)
rand.eff=formula(paste("~1|",sub.ID))

if (is.na(ped$theta)) { fit <- try(lmekin(fixed=fix.eff,data=x,random = rand.eff,varlist=list(kmat))) } else { 
if ( ped$theta<=0 ) fit <- try(mlekin.get.theta(fixed=fix.eff,data=x,random = rand.eff,varlist=list(kmat)))
if ( ped$theta >0 ) fit <- try(mlekin.fixed.theta(fixed=fix.eff,data=x,random = rand.eff,varlist=list(kmat),fixed.theta=theta))
}

if (class(fit)!="try-error")
{

if (!is.na(ped$theta)) { if  (ped$theta<=0)   fit$optim.theta -> theta } 

if (xi %in% rownames(fit$ctable)) 
{
fit$ctable[xi,1]->Beta
fit$ctable[xi,2]->SE
fit$ctable[xi,3]->t
fit$ctable[xi,4]->p.value
sampleN=fit$n
DF=fit$df.res
}
Pval0=cor.test(x[,xi],x[,yi])$p.value
cor(x[,xi],x[,yi],use="pairwise.complete.obs")^2 -> rsquare
Beta=signif(Beta,4)
SE=signif(SE,4)
t=signif(t,4)
p.value=signif(p.value,4)
Pval0=signif(Pval0,4)
rsquare=signif(rsquare,4)
} # if no fit error
} # nxy >10
pp=rbind(pp,cbind(markname,phenotype, as.data.frame(cbind(Beta,SE,t,p.value,DF,sampleN,Pval0,rsquare,theta)) ))
}} #yi,xi
pp
} # end kin.test


############ QQ plot #################################

myqqplot=function(ps,ylm=NULL)
{
as.data.frame(ps)->ps
for (pname in colnames(ps))
{
pi=ps[,pname];pi=pi[!is.na(pi)];pi=-log10(pi);pi=sort(pi)
p0=c(1:length(pi))/length(pi);p0=-log10(p0);p0=sort(p0)
pp=cbind(p0,pi)
matplot(p0,pp,main=pname,ylab="",xlab="",ylim=ylm,type="l",lty=1,lwd=2)
#legend(x="bottomright",legend=colnames(x),col=co,lty=1,lwd=2,box.lty=0)
}
}

############## ROC ############################
myroc=function(tt,ps,effcol,fpr=NULL)
{
if (is.null(fpr)) fpr=c((1:9)/100,(1:10)/10)
ttp=as.data.frame(fpr)
for (gi in ps)
{
p=tt[,gi]
p[is.na(p)]=1
p0=p[tt[,effcol]==0]
qk=quantile(p0,fpr)
pi=p[tt[,effcol]!=0]
ni=length(pi)
powi=NULL;for (j in 1:nrow(ttp)) powi=c(powi,signif(sum(pi<=qk[j])/ni,5))
ttp=cbind(ttp,powi)
}
colnames(ttp)=c("FPR",ps)
ttp
}

####### ROC PLOT ################

mymatplot=function(x,y=NULL,xlm=c(0,1),ylm=c(0,1),lgs=NULL,tit="",cos=c("black","blue","red","green","brown","pink","purple"))
{
if (is.null(y)) {y=x[,-1];colnames(y)=colnames(x)[-1];x=x[,1]}
if (is.null(lgs)) lgs=colnames(y)
matplot(x,y,xlim=xlm,ylim=ylm,ylab="TPR",xlab="FPR",main=tit,cex.main=2,
type="l",lwd=2,lty=1,pch="*",cex.lab=1.5,cex.axis=1.5,col=cos)
legend(x="bottomright",legend=lgs,lwd=2,lty=1,col=cos,cex=1.5,bty="n")
}

####################### get sum weights Madson & Browning ################################

sumweights=function (x,ctrid=NULL)
{
#x: genotyps 0,1,2; 1,2 are rare/minor
#ctrid: control id

xii=x
if (is.null(ctrid)) ctrid=c(1:nrow(xii))
cid=rep(F,nrow(xii))
cid[ctrid]=T

###rare allele total counts
x0=rowSums((xii==2),na.rm=T)
x1=rowSums((xii==1),na.rm=T)
x2=x0+x1
xt=2*x0+x1
###Collaping
z0=as.numeric(x0>0)
z1=as.numeric(x1>0)
z2=as.numeric(x2>0)

###weighted sum
sum(xt[cid],na.rm=T)-> xc1
sum(xt,na.rm=T)-xc1 -> xc2
sum(z2[cid],na.rm=T)-> nc1
sum(z2,na.rm=T)-nc1 -> nc2

cid[is.na(cid)]=F
niu=sum(cid) # control number
miu=colSums(xii[cid,],na.rm=T) # rare var number per snp in controls
qi=(miu+1)/(2*niu+2)
wi=sqrt(nrow(xii)*qi*(1-qi))
#Iij=xii*NA;Iij[xii==0]=1;Iij[xii==1]=0;Iij[xii==2]=0
#rowSums(t(t(Iij)/wi),na.rm=T)->w0
wi=1/wi
tt=list(wi=wi,xc1=xc1,xc2=xc2,nc1=nc1,nc2=nc2)
tt
}


##########################################
g.simu.unif=function(n,m,f,name)
{
    g.prob=c((1-f)^2,2*f*(1-f),f^2)
    sample(x=c(0,1,2),size=m*n,replace=T,prob=g.prob)->temp
    matrix(temp,n,m)->temp
    colnames(temp)=paste(name,1:m,sep="")
    temp
}

##########################################
g.simu.varf=function(n,m,f,fsd,name)
{
    fs=abs(rnorm(m,f,fsd))
    tt=NULL
    for (fi in fs)
    {
        g.prob=c((1-fi)^2,2*fi*(1-fi),fi^2)
        sample(x=c(0,1,2),size=n,replace=T,prob=g.prob)->temp
        tt=cbind(tt,temp)
    }
    colnames(tt)=paste(name,1:m,sep="")
    tt
}


### genotype(x) simulation ###
genotype.simu=function(n.pool,n.mk,r.up,r.dw,f.ra,f.sd=0,vtype="number")
{
    N=n.pool

    if (vtype=="percent")
    {
        x.up=NULL; M=as.integer(r.up*n.mk+0.5); if (M>0) g.simu(n=N,m=M,f=f.ra,name="up")->x.up
        x.dw=NULL; M=as.integer(r.dw*n.mk+0.5); if (M>0) g.simu(n=N,m=M,f=f.ra,name="dw")->x.dw
        x.md=NULL; M=n.mk-as.integer(r.up*n.mk+0.5)-as.integer(r.dw*n.mk+0.5)
        if (M>0) g.simu(n=N,m=M,f=f.ra,name="md")->x.md
    }

    if (vtype=="number")
    {
        x.up=NULL; M=r.up; if (M>0) g.simu.varf(n=N,m=M,f=f.ra,fsd=f.sd,name="up")->x.up
        x.dw=NULL; M=r.dw; if (M>0) g.simu.varf(n=N,m=M,f=f.ra,fsd=f.sd,name="dw")->x.dw
        x.md=NULL; M=n.mk-r.up-r.dw
        if (M>0) g.simu.varf(n=N,m=M,f=f.ra,fsd=f.sd,name="md")->x.md
    }

    x=cbind(x.up,x.dw,x.md)
    rownames(x)=c(1:nrow(x))
    x
}

###########################################################
phenotype.simu=function(x,eff.ra,eff.sd,model="additive")
{
if (model=="additive")
{
eff=rnorm(ncol(x),eff.ra,eff.sd)
grep("up",colnames(x))-> col.up
grep("dw",colnames(x))-> col.dw
grep("md",colnames(x))-> col.md
if (length(col.up)>0) eff[col.up]=abs(eff[col.up])
if (length(col.dw)>0) eff[col.dw]=-abs(eff[col.dw])
if (length(col.md)>0) eff[col.md]=0
y=rnorm(nrow(x)) +  colSums(t(x)*eff)
}
y
}

################################################
sampling.simu=function(y0,n.sample,design="RDM")
{
    if (design=="RDM")
    {
        id=sample(1:length(y0),n.sample)
        y=y0[id]
        grp=NULL
    }

    if (design=="LR")
    {
        id=1:length(y0)
        n.sample/length(y0)->pct
        id[y<=quantile(y,pct/2)] -> Lid
        id[y>=quantile(y,(1-pct/2))] -> Rid
        id=c(Lid,Rid)
        y=y0[id]
        grp=c(rep("L",length(Lid)),rep("R",length(Rid)))
    }

    if (design=="LCR")
    {
        id=1:length(y0)
        n.sample/length(y0)->pct
        id[y<=quantile(y,pct/3)] -> Lid
        id[y>=quantile(y,(1-pct/3))] -> Rid
        id[y>=quantile(y,(0.5-pct/3/2)) & y<=quantile(y,(0.5+pct/3/2))] -> Cid
        id=c(Lid,Cid,Rid)
        y=y0[id]
        grp=c(rep("L",length(Lid)),rep("C",length(Cid)),rep("R",length(Rid)))
    }

    if (design=="LC")
    {
        id=1:length(y0)
        n.sample/length(y0)->pct
        id[y<=quantile(y,pct/2)] -> Lid
        id[y>=quantile(y,(0.5-pct/2/2)) & y<=quantile(y,(0.5+pct/2/2))] -> Cid
        id=c(Lid,Cid)
        y=y0[id]
        grp=c(rep("L",length(Lid)),rep("C",length(Cid)))
    }

    if (design=="CR")
    {
        id=1:length(y0)
        n.sample/length(y0)->pct
        id[y>=quantile(y,(1-pct/2))] -> Rid
        id[y>=quantile(y,(0.5-pct/2/2)) & y<=quantile(y,(0.5+pct/2/2))] -> Cid
        id=c(Cid,Rid)
        y=y0[id]
        grp=c(rep("C",length(Cid)),rep("R",length(Rid)))
    }

    samp=list(y=y,id=id,grp=grp,design=design)
    samp
}

######################## variable selection

vSelect=function(x,y,w=NULL,method="hoffmann")
{
x=as.data.frame(x)

if (ncol(x)==0) z=list(id=0,chi=NA)
if (ncol(x)==1) z=list(id=1,chi=NA)
if (ncol(x)>1) 
{
if (method=="hoffmann") #step-up,by Hoffmann,PLOSone,2010,5(11)p3
{

###
#z=samp
#z$x=z$x[,sd(z$x,na.rm=T)>0]
#x=z$x
#y=z$y
#x=as.data.frame(x)
if (is.null(w))
{
w=as.numeric(cor(x,y))
w[w>0]=1
w[w<0]=-1
}
###

xkm=colMeans(x,na.rm=T)
ym=mean(y,na.rm=T)
t((t(x)-xkm)*w)*(y-ym)->u
ks0=1:ncol(x)
kid=NULL
chi00=0
tt=NULL
while(length(kid)<length(ks0))  #----------------
{

if (is.null(kid)) ks=ks0
if (!is.null(kid)) ks=ks0[-kid]
chi0=0
ki=0
for (k in ks)
{ 
if (sd(u[,k],na.rm=T)>0)
{
c(kid,k)->ksi
sum(u[,ksi],na.rm=T)^2 / sum(rowSums(cbind(u[,ksi],0),na.rm=T)^2,na.rm=T) -> chik
if (chik>chi0) {chi0=chik;ki=k}
}
}

if (chi0>=chi00) {kid=c(kid,ki);chi00=chi0}
#kid=c(kid,ki)
#tt=c(tt,chi0)
#print(kid)
#print(chi0)
if (chi0<chi00) break

} # while -------------------------------------------

ks0*0->sele
sele[kid]=1
if (chi00==0) chi00=NA
z=list(id=sele,chi=chi00)

} # hoffman

} # ncol>1

z

}

###########################################



############# kinship model to estimate optim.theta ################

mlekin.get.theta= function (fixed, data = parent.frame(), random, varlist = NULL, 
    variance, sparse = c(20, 0.05), rescale = T, pdcheck = T, 
    subset, weight, na.action) 
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    if (missing(variance)) 
        theta <- NULL
    else theta <- variance
    reSt <- reStruct(random, REML = F, data = NULL)
    gform <- getGroupsFormula(reSt)
    if (is.null(gform)) {
        temp.fixed <- fixed
        gvars <- NULL
    }
    else {
        gvars <- all.vars(random)
        fvars <- all.vars(formula)
        gvars <- gvars[is.na(match(gvars, fvars))]
        temp.fixed <- paste(deparse(as.vector(fixed)), collapse = "")
        temp.fixed <- paste(temp.fixed, paste(gvars, collapse = "+"), 
            sep = "+")
        temp.fixed <- as.formula(temp.fixed)
    }
    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- terms(fixed)
    X <- model.matrix(Terms, m)
    Y <- model.extract(m, "response")
    n <- length(Y)
    weights <- model.extract(m, "weights")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, n)
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    ncluster <- length(gvars)
    if (ncluster == 0) 
        stop("No grouping variables found")
    groups <- getGroups(m, gform)
    temp <- coxme.varcheck(ncluster, varlist, n, gvars, groups, 
        sparse, rescale, pdcheck)
    varlist <- temp$varlist
    kindex <- temp$kindex
    ntheta <- temp$ntheta
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i] == 1) 
            theta.names <- c(theta.names, gvars[i])
        else theta.names <- c(theta.names, paste(gvars[i], 1:ntheta[i], 
            sep = ""))
    }
    if (length(theta) == 0) 
        theta <- rep(0, sum(ntheta))
    else if (length(theta) != sum(ntheta)) 
        stop("Wrong length for theta")
    names(theta) <- theta.names
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i] == 1) 
            theta.names <- c(theta.names, gvars[i])
        else theta.names <- c(theta.names, paste(gvars[i], 1:ntheta[i], 
            sep = ""))
    }
    if (length(theta) == 0) 
        theta <- rep(0, sum(ntheta))
    else if (length(theta) != sum(ntheta)) 
        stop("Wrong length for theta")
    names(theta) <- theta.names
    tindex <- which(theta == 0)
    if (ncluster > 1) 
        stop("function can have only 1 random effect")
    varlist <- varlist[[1]]
    kindex <- kindex[, 1]
    if (max(kindex) != n) 
        stop("The random effect must be 1 per subject")
    ntheta <- ntheta[1]
    kindex2 <- integer(n)
    kindex2[kindex] <- 1:n
    logfun <- function(itheta, X, Y, varlist, theta, tindex, 
        center) {
        theta[tindex] <- exp(itheta)
        tkmat <- varlist[[1]]
        tkmat@blocks <- tkmat@blocks * theta[1]
        diag(tkmat) <- diag(tkmat + 1)
        if (length(varlist) > 1) {
            for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks * 
                theta[i] + tkmat@blocks
        }
        tkmat@blocks <- tkmat@blocks/tkmat@blocks[1]
        gk <- gchol(tkmat)
        newx <- solve(gk, X, full = FALSE)
        newy <- solve(gk, Y, full = FALSE)
        resid <- qr.resid(qr(newx), newy)
        n <- length(Y)
        loglik <- (n/2) * (log(mean(resid^2)) - center) + sum(log(diag(gk)))/2
        loglik
    }
    newX <- X[kindex2, ]
    newY <- as.vector(Y[kindex2])
    dimnames(newX) <- NULL
    if (length(tindex) > 0) {
        center <- log(mean((Y - mean(Y))^2))
        nfit <- optim(par = rep(-1, length(tindex)), logfun, 
            method = "L-BFGS-B", lower = log(1e-05), X = newX, 
            Y = newY, varlist = varlist, theta = theta, tindex = tindex, 
            center = center)
        iter <- nfit$counts
        theta[tindex] <- exp(nfit$par) ;  optim.theta=theta 
    }
    else iter <- 0
    tkmat <- varlist[[1]]
    tkmat@blocks <- tkmat@blocks * theta[1]
    diag(tkmat) <- diag(tkmat + 1)
    if (length(varlist) > 1) {
        for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks * 
            theta[i] + tkmat@blocks
    }
    gk <- gchol(tkmat)
    xok <- as.matrix(solve(gk, newX, full = F))
    yok <- solve(gk, newY, full = FALSE)
    lfit <- lm(yok ~ 0 + xok)
    names(lfit$coefficients) <- dimnames(X)[[2]]
    ls <- summary(lfit)
    resid.var <- mean(lfit$residuals^2)
    theta <- c(theta * resid.var, resid.var)
    names(theta) <- c(theta.names, "resid")
    fitted <- c(X %*% lfit$coef)
    residuals <- Y - fitted
    frail <- residuals[kindex2]
    names(frail) <- groups
    fcoef <- lfit$coef
    call$fixed <- fixed
    call$random <- random
    fit <- list(coefficients = list(fixed = fcoef, random = frail), 
        theta = theta, variance = ls$cov.unscaled * ls$sigma^2, 
        ctable = ls$coefficients, residuals = residuals, fitted.values = fitted, 
        effects = lfit$effects, rank = lfit$rank, assign = lfit$assign, 
        df.residual = lfit$df.residual - length(theta), loglik = (-n/2) * 
            (log(mean(lfit$residuals^2)) + 1 + log(2 * pi)) - 
            sum(log(diag(gk)))/2, iter = iter, n = n, call = call, 
        method = "ML", optim.theta=optim.theta)
    na.action <- attr(m, "na.action")
    if (length(na.action)) 
        fit$na.action <- na.action
    oldClass(fit) <- c("lmekin")
    fit
}

################ kinship model with known theta #############


mlekin.fixed.theta= function (fixed, data = parent.frame(), random, varlist = NULL, 
    variance, sparse = c(20, 0.05), rescale = T, pdcheck = T, 
    subset, weight, na.action, fixed.theta) 
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    if (missing(variance)) 
        theta <- NULL
    else theta <- variance
    reSt <- reStruct(random, REML = F, data = NULL)
    gform <- getGroupsFormula(reSt)
    if (is.null(gform)) {
        temp.fixed <- fixed
        gvars <- NULL
    }
    else {
        gvars <- all.vars(random)
        fvars <- all.vars(formula)
        gvars <- gvars[is.na(match(gvars, fvars))]
        temp.fixed <- paste(deparse(as.vector(fixed)), collapse = "")
        temp.fixed <- paste(temp.fixed, paste(gvars, collapse = "+"), 
            sep = "+")
        temp.fixed <- as.formula(temp.fixed)
    }
    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- terms(fixed)
    X <- model.matrix(Terms, m)
    Y <- model.extract(m, "response")
    n <- length(Y)
    weights <- model.extract(m, "weights")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, n)
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    ncluster <- length(gvars)
    if (ncluster == 0) 
        stop("No grouping variables found")
    groups <- getGroups(m, gform)
    temp <- coxme.varcheck(ncluster, varlist, n, gvars, groups, 
        sparse, rescale, pdcheck)
    varlist <- temp$varlist
    kindex <- temp$kindex
    ntheta <- temp$ntheta
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i] == 1) 
            theta.names <- c(theta.names, gvars[i])
        else theta.names <- c(theta.names, paste(gvars[i], 1:ntheta[i], 
            sep = ""))
    }
    if (length(theta) == 0) 
        theta <- rep(0, sum(ntheta))
    else if (length(theta) != sum(ntheta)) 
        stop("Wrong length for theta")
    names(theta) <- theta.names
    theta.names <- NULL
    for (i in 1:ncluster) {
        if (ntheta[i] == 1) 
            theta.names <- c(theta.names, gvars[i])
        else theta.names <- c(theta.names, paste(gvars[i], 1:ntheta[i], 
            sep = ""))
    }
    if (length(theta) == 0) 
        theta <- rep(0, sum(ntheta))
    else if (length(theta) != sum(ntheta)) 
        stop("Wrong length for theta")
    names(theta) <- theta.names
    tindex <- which(theta == 0)
    if (ncluster > 1) 
        stop("function can have only 1 random effect")
    varlist <- varlist[[1]]
    kindex <- kindex[, 1]
    if (max(kindex) != n) 
        stop("The random effect must be 1 per subject")
    ntheta <- ntheta[1]
    kindex2 <- integer(n)
    kindex2[kindex] <- 1:n
    logfun <- function(itheta, X, Y, varlist, theta, tindex, 
        center) {
        theta[tindex] <- exp(itheta)
        tkmat <- varlist[[1]]
        tkmat@blocks <- tkmat@blocks * theta[1]
        diag(tkmat) <- diag(tkmat + 1)
        if (length(varlist) > 1) {
            for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks * 
                theta[i] + tkmat@blocks
        }
        tkmat@blocks <- tkmat@blocks/tkmat@blocks[1]
        gk <- gchol(tkmat)
        newx <- solve(gk, X, full = FALSE)
        newy <- solve(gk, Y, full = FALSE)
        resid <- qr.resid(qr(newx), newy)
        n <- length(Y)
        loglik <- (n/2) * (log(mean(resid^2)) - center) + sum(log(diag(gk)))/2
        loglik
    }
    newX <- X[kindex2, ]
    newY <- as.vector(Y[kindex2])
    dimnames(newX) <- NULL
    if (length(tindex) > 0) {
        center <- log(mean((Y - mean(Y))^2))
        #nfit <- optim(par = rep(-1, length(tindex)), logfun, 
        #    method = "L-BFGS-B", lower = log(1e-05), X = newX, 
        #    Y = newY, varlist = varlist, theta = theta, tindex = tindex, 
        #    center = center)
        #iter <- nfit$counts
        #theta[tindex] <- exp(nfit$par) ;  optim.theta=theta 
    
iter=0;theta[tindex]=fixed.theta; # optim.theta=theta 

}
    else iter <- 0
    tkmat <- varlist[[1]]
    tkmat@blocks <- tkmat@blocks * theta[1]
    diag(tkmat) <- diag(tkmat + 1)
    if (length(varlist) > 1) {
        for (i in 2:length(varlist)) tkmat@blocks <- varlist[[i]]@blocks * 
            theta[i] + tkmat@blocks
    }
    gk <- gchol(tkmat)
    xok <- as.matrix(solve(gk, newX, full = F))
    yok <- solve(gk, newY, full = FALSE)
    lfit <- lm(yok ~ 0 + xok)
    names(lfit$coefficients) <- dimnames(X)[[2]]
    ls <- summary(lfit)
    resid.var <- mean(lfit$residuals^2)
    theta <- c(theta * resid.var, resid.var)
    names(theta) <- c(theta.names, "resid")
    fitted <- c(X %*% lfit$coef)
    residuals <- Y - fitted
    frail <- residuals[kindex2]
    names(frail) <- groups
    fcoef <- lfit$coef
    call$fixed <- fixed
    call$random <- random
    fit <- list(coefficients = list(fixed = fcoef, random = frail), 
        theta = theta, variance = ls$cov.unscaled * ls$sigma^2, 
        ctable = ls$coefficients, residuals = residuals, fitted.values = fitted, 
        effects = lfit$effects, rank = lfit$rank, assign = lfit$assign, 
        df.residual = lfit$df.residual - length(theta), loglik = (-n/2) * 
            (log(mean(lfit$residuals^2)) + 1 + log(2 * pi)) - 
            sum(log(diag(gk)))/2, iter = iter, n = n, call = call, 
        method = "ML", fixed.theta=fixed.theta)
    na.action <- attr(m, "na.action")
    if (length(na.action)) 
        fit$na.action <- na.action
    oldClass(fit) <- c("lmekin")
    fit
}


########################################## maf()
# Genotype : MAF
maf=function(x)
{
colSums(x==2,na.rm=T)->aa
colSums(x==1,na.rm=T)->ab
colSums(x==0,na.rm=T)->bb
f=(aa+ab/2)/(aa+ab+bb)
f
}
#########

dqrls_wrapper <- function(X, y) {
    rnames <- colnames(X);

    n <- dim(X)[1];
    p <- dim(X)[2];

    fit <- .Fortran("dqrls",
                    qr = X, n = n,
                    p = p, y = y, ny = 1L,
                    tol = 1e-7,
                    coefficients = double(p),
                    residuals = double(n),
                    effects = double(n),
                    rank = integer(1L),
                    pivot = 1L:p, qraux = double(p),
                    work = double(2 * p),
                    PACKAGE = "base");

    # degenerate matrix! some columns were pivoted out, we'll need to ignore
    # them.
    if (fit$rank < p) {
        X <- X[,fit$pivot][,1:fit$rank, drop=FALSE];
        rnames <- colnames(X);
    }

    # grab the coefficients, X %*% theta is our linear predictor
    theta <- fit$coeff[1:fit$rank];

    # compute (X^t X)^{-1}; approx covariance matrix
    XtXinv <- chol2inv(fit$qr[1:fit$rank,1:fit$rank,drop=FALSE]);

    # let's make a p-value
    df.residual <- dim(X)[1] - fit$rank;
    dispersion <- sum( (X %*% theta - y)^2 ) / df.residual;
    s.err <- sqrt(dispersion * diag(XtXinv));
    tvalue <- theta / s.err;
    pvalue <- 2 * pt(-abs(tvalue), df.residual);

    # build return value
    coeff <- cbind(theta, s.err, tvalue, pvalue);
    dimnames(coeff) <- list(
        rnames,
        c("Estimate", "Std.Error", "t value", "Pr(>|t|)"));
    fit <- list(coeffs = coeff, XtXinv = XtXinv);

    fit;
}

#################### myglm
myglm_orig=function(z,trait,variant,covar=NULL,ytype)
{
    if (is.null(covar)) {
        model=formula(paste(trait,"~",variant)) 
    } else {
        model=formula(paste(trait,"~",variant,"+",covar))
    }
    if (ytype=="B") fit=try(summary(glm(formula=model,data=z,family=binomial(link = "logit"))))
    if (ytype=="Q") fit=try(summary(lm(formula=model,data=z,family=gaussian(link = "identity"))))
    fit$p.value=fit$coef[variant,4]
    fit
}
###################

myglm=function(z,trait,variant,covar=NULL,ytype,covar_matrix=NULL)
{
    model <- NULL

    # If a precomputed covariant matrix was passed in, use that, otherwise,
    # generate one (if there are in fact covariates).
    if (ytype == "Q") { 
        if (is.null(covar_matrix) && !is.null(covar)) {
            covar_matrix <- model.frame(formula=formula(paste(trait,"~",covar,"-1")), z);
            mt <- attr(covar_matrix, "terms") # allow model.frame to have updated it
            covar_matrix <- model.matrix(mt, covar_matrix, NULL);
        }
        model_txt = paste(trait,"~",variant)
        model=formula(model_txt)
    } else {
        if (is.null(covar)) {
            model=formula(paste(trait,"~",variant))
        } else {
            model=formula(paste(trait,"~",variant,"+",covar))
        }
    }


    fit <- NULL;
    if (ytype=="B") fit=try(summary(glm(formula=model,data=z,family=binomial(link = "logit"))))
    if (ytype=="Q") {

        mf <- model.frame(formula=model, data=z);
        X <- model.matrix(model, mf, NULL);
        y <- model.response(mf, "any") # factors are allowed as response vars
        fit = NULL
        X <- cbind(X, covar_matrix);
        fit = dqrls_wrapper(X, y);
    }

    fit$p.value=fit$coef[variant,4]
    fit
}
