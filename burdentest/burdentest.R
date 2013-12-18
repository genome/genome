#!/usr/bin/env Rscript

rm(list=ls())
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(script.name)

params = commandArgs(trailingOnly=TRUE)

option.file=gsub(" ","",params[1])  #option.file="option_file_asms"
print(paste("Sourcing option file", option.file))
source(option.file)

# set up multi-core processing
library(doMC);
registerDoMC(cores=num.cores);


library.file=paste(script.dir, "rarelib_skat.R", sep="/");
#option.file="/gscuser/qzhang/gstat/burdentest/option_file_asms_burden"
source(library.file)


trait.type=gsub(" ","",params[2]) #trait.type="Q"
trait=gsub(" ","",params[3])      #trait="trigRES"
gene=gsub(" ","",params[4])      #gene="ABCA1"

permu.par=gsub(" ","",params[5])      #permu.par="10000:1:1"

permu.par=strsplit(permu.par,split=":")[[1]]
permu.n=as.numeric(permu.par[1])
rd.seed=as.numeric(permu.par[2])
#if (is.na(rd.seed)) rd.seed=0
if (rd.seed>0) set.seed(rd.seed)
pp.n=as.numeric(permu.par[3])
#if (is.na(pp.n)) pp.n=1

add.cov=gsub(" ","",params[6])
if (is.na(add.cov)) add.cov="NONE"
cova=covariates
if (cova!="NONE" & add.cov!="NONE")  covariates=paste(covariates,add.cov,sep="+")
if (cova!="NONE" & add.cov=="NONE")  covariates=cova
if (cova=="NONE" & add.cov!="NONE")  covariates=add.cov
if (cova=="NONE" & add.cov=="NONE")  covariates=NULL

print(library.file)
print(option.file)
print(trait.type)
print(trait)
print(gene)
print(permu.n)
print(rd.seed)
print(pp.n)
print(covariates)
print(genotype.file)
print(phenotype.file)
print(anno.file)


read.table(genotype.file,head=T,sep=gfile.delimiter,na.strings=c("NA",".","",missing.data))->x
read.table(phenotype.file,head=T,sep=pfile.delimiter,na.strings=c("NA",".","",missing.data))->y

lstx=NA;lsty=NA
if (samplelist.dir!="NONE")
{
read.table(paste(samplelist.dir, paste(trait,"snp",sep=".") ,sep="/"),header=F)->lstx
read.table(paste(samplelist.dir, paste(trait,"people",sep=".") ,sep="/"),header=T)->lsty
lstx=as.character(lstx[,1])
names(lsty)=toupper(names(lsty))
lsty=as.character(lsty[lsty[,toupper(gene)]==1,1])
}


m <- read.table(
    anno.file,
    check.names=F,
    comment.char="",
    head=T,
    sep=gfile.delimiter,
    na.strings=c("NA",
    ".",
    "",
    missing.data))
# strip optional leading # from header
colnames(m)[1] <- gsub("^#", "", colnames(m)[1])

if (vtype.use[1]!="ALL") m=m[(m[,vtype.col] %in% vtype.use),]

############################## Batch least squares class 
bls.init <- function(data, variantNames=NULL) {
    if (is.null(variantNames)) {
        variantNames = colnames(data)
    }

    bls <- list(
        variantNames = variantNames,
        cache = list()
        )
    class(bls) <- "bls"
    numSamples <- dim(data)[1]

    for (name in variantNames) {
        X <- cbind(rep(1, numSamples), data[,name])
        qrd <- qr(X)
        covarianceBase <- chol2inv(qrd$qr[1:qrd$rank,1:qrd$rank,drop=FALSE]);
        leftInverse <- covarianceBase %*% t(X)
        bls$cache[[name]] = list(rank=qrd$rank, covarianceBase=covarianceBase, leftInverse=leftInverse, X=X)
    }
    bls
}

bls.pvalues <- function(bls, y, variantNames = NULL) {
    if (class(bls) != "bls") {
        stop(paste("Object passed to bls.pvalues is not of class bls:",class(bls)))
    }
    if (is.null(variantNames)) {
        variantNames = bls$variantNames
    }

    result = mat.or.vec(length(bls$variantNames),  4)
    rownames(result) = bls$variantNames
    colnames(result) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)")

    i = 1
    for (name in variantNames) {
        cacheEntry = bls$cache[[name]]
        if (is.null(cacheEntry)) {
            stop(paste("Failed to find cached bls cache entryfor", name))
        }
        covarianceBase = cacheEntry$covarianceBase
        leftInv = cacheEntry$leftInverse
        X = cacheEntry$X
        rank = cacheEntry$rank

        theta <- leftInv %*% y
        # let's make a p-value
        df.residual <- dim(X)[1] - rank;
        residuals = X %*% theta - y
        dispersion <- sum( residuals^2 ) / df.residual;
        s.err <- sqrt(dispersion * diag(covarianceBase));
        tvalue <- theta / s.err;
        pvalue <- 2 * pt(-abs(tvalue), df.residual);

        # build return value
        # we are interested in the 2nd row only, the first is the bias term
        row = cbind(theta, s.err, tvalue, pvalue)[2,]
        result[i,] = row
        i = i + 1
    }
    result
}

######################################### collapsing.test()

collapsing.test=function(
    x,y,m,trait,gene,
    ytype,covar=NULL,class_unit,gvid,avid,gsubid,psubid,maf_cutoff=0.01,permun=10000,ped=NULL,lstx=NA,lsty=NA
    )
{

    ### select data for analysis
    xs=as.character(m[m[,class_unit]==gene,avid])
    if (!is.na(lstx)) xs=xs[xs %in% lstx]

    if (!is.null(covar)) {covars=strsplit(covar,split="[+]")[[1]]; covars=unique(covars)}
    if (is.null(covar)) covars=NULL

    rst=NULL
    if (length(xs)>1)
    {

        if (gsubid=="FIRST_ROW")
        {
            x[x[,gvid] %in% c(xs),] -> x
            rownames(x)=x[,gvid]
            x[,colnames(x)!=gvid]->x
            x=as.data.frame(t(x))
            x$FIRST_ROW=rownames(x)
            y[,psubid] = gsub("-",".",y[,psubid])
        }

        allid=y[,psubid] # for fixed permutation retaining correlation between traits

        colnames(x) %in% xs->xsid
        colnames(x)[xsid]=paste("V",colnames(x)[xsid],sep="")
        xs=colnames(x)[xsid]

        y[,colnames(y) %in% c(psubid,trait,covars)]->y
        y=y[rowSums(is.na(y))==0,]
        if (!is.na(lsty)) y=y[y[,psubid] %in% gsub("-",".",lsty),]

        x[,colnames(x) %in% c(gsubid,xs)]->x

        xy=merge(x,y,by.x=gsubid,by.y=psubid)


        ### build data object for analysis
        design="RDM"
        id=xy[,gsubid]
        grp=NA
        yi=xy[,trait]
        xi=xy[,colnames(xy) %in% xs]

        ###  MAF
        maf=function(x)
        {
            colSums(x==2,na.rm=T)->aa
            colSums(x==1,na.rm=T)->ab
            colSums(x==0,na.rm=T)->bb
            f=(aa+ab/2)/(aa+ab+bb)
            f
        }
        ###

        f=maf(xi)
        if (sum(f>0.5,na.rm=T)>0) { xi[,f>0.5]=2-xi[,f>0.5]; f[f>0.5]=1-f[f>0.5] }
        sele = (f>0 & f<=maf_cutoff & !is.na(f))

        if (sum(sele)>1)
        {
            xi=xi[,sele];names(xi)->xs
            xy[,xs]=xi
            xy=xy[rowSums(is.na(xy[,xs]))==0,]
            yi=xy[,trait]
            xi=xy[,xs]
            id=xy[,gsubid]

            xs[sd(xi,na.rm=T)>0]->xs

            if (length(xs)>1)
            {

if (permun<0)
{
permun=abs(permun)
sample(1:nrow(xy))->outer.permu.id
xy[,xs]=xy[outer.permu.id,xs]
}

                # batch least squares
                bls = NULL
                if (!is.null(covars) && ytype == 'Q') {
                    print("Pre-fitting covariate model")
                    print(head(xy[,c(trait,covars)]))
                    model=formula(paste(trait,"~",covar)) 
                    fit <- lm(formula=model,data=xy[,c(trait,covars)])
                    print("Covar model")
                    print(summary(fit))
                    resid <- fit$resi
                    xy[,trait] = resid
                    covar = covars = NULL
                    bls <- bls.init(xy[, xs], xs)
                }


                z=list(
                    design=design,id=id,allid=allid,grp=grp,ped=ped,ytype=ytype,data=xy[,c(trait,covars,xs)],
                    trait=trait,variants=xs,covar=covar,covars=covars,bls=bls
                )
                collapse.test(z,permun)->rst
                rst$gene=gene
                rst$maf=maf_cutoff
                rst$V=length(xs)
                rst$N=nrow(xy)

            } # final xs>1
        } #sele>1

    } #(length(xs)>1)

    rst

}

##################################################
ppp=NULL

tm_start <- proc.time();
ppp=NULL
for (ppi in c(1:pp.n))
{
    print(paste("p-value permutation", ppi));

    try(collapsing.test(x=x,y=y,m=m,
    trait=trait,gene=gene,ytype=trait.type,
    covar=covariates,
    gvid=gfile.vid,
    avid=afile.vid,
    class_unit=gene.col,
    gsubid=gfile.sid,
    psubid=pfile.sid,
    maf_cutoff=maf.cutoff,
    permun=permu.n,
    ped=NULL,
    lstx=lstx,
    lsty=lsty))->rst

    if (!file.exists(out.dir)==T) dir.create(out.dir,recursive = FALSE)

    if (class(rst)=="try-error") writeLines(rst,con=paste(out.dir,"/",trait,"_",gene,".error",sep=""))
    if (is.null(rst)) write.table(rst,file=paste(out.dir,"/",trait,"_",gene,".null",sep=""))
    if (!is.null(rst))
    {

        rst$z$trait->Trait
        rst$gene->Gene
        rst$maf->MAF
        rst$V->V
        rst$N->N
        CMC=rst$p["C01"]
        pCMC=rst$pp["C01"]
        WSS=rst$p["WSS"]
        pWSS=rst$pp["WSS"]
        aSum=rst$pp["aSum"]
        PWST=rst$pp["PWST"]
        SPWST=rst$pp["SPWST"]
        SPWST.up=rst$pp["xuplog"]
        SPWST.down=rst$pp["xdwlog"]

        SKAT=rst$p["skat"]
        pSKAT=rst$pp["skat"]
        MREG=rst$p["mreg"]
        pMREG=rst$pp["mreg"]

        vt.maf=rst$p["VT_MAF"]
        VT=rst$pp["VT"]

        pp=cbind(Trait,Gene,N,V,MAF,CMC,pCMC,WSS,pWSS,aSum,PWST,SPWST,SPWST.up,SPWST.down,SKAT,pSKAT,MREG,pMREG,VT,vt.maf)
        #ppp=rbind(ppp,cbind(ppi,pp))

        if (ppi==1)
        {
            
	     #burden test output
            write.csv(pp,file=paste(out.dir,"/",trait,"_",gene,".burden.csv",sep=""),quote=F,row.names=F)
            
            #single-var test output
            p=NA
            p=rst$svt
            Variant=rownames(p)
            colnames(p)[1:4]=c("Beta","SE","t","P")
            p=signif(p,5)
            p=cbind(Variant,p)
            write.csv(p,file=paste(out.dir,"/",trait,"_",gene,".single.csv",sep=""),quote=F,row.names=F)
        }

    } # rst not null

} # ppi loop

tm_end <- proc.time();
print("Running time:");
print(tm_end - tm_start);

#write.csv(ppp,file=paste(out.dir,"/",trait,"_",gene,".ppi.csv",sep=""),quote=F,row.names=F)

############################################################################################

temp4debug=function(x)
{

x=x
y=y
m=m
trait=trait
gene=gene
ytype=trait.type
covar=covariates
gvid=gfile.vid
avid=afile.vid
class_unit=gene.col
gsubid=gfile.sid
psubid=pfile.sid
maf_cutoff=maf.cutoff
permun=permu.n
ped=NULL
lstx=lstx
lsty=lsty

}


