rm(list=ls())

library.file="/gscuser/qzhang/gstat/burdentest/rarelib.R"
#option.file="/gscuser/qzhang/gstat/burdentest/option_file_asms_burden"
source(library.file)

option.file=gsub(" ","",commandArgs()[3])
source(option.file)

trait.type=gsub(" ","",commandArgs()[4])
trait=gsub(" ","",commandArgs()[5])
gene=gsub(" ","",commandArgs()[6])
permu.n=gsub(" ","",commandArgs()[7]); permu.n=as.numeric(permu.n)

print(library.file)
print(option.file)
print(trait.type)
print(trait)
print(gene)
print(permu.n)

read.table(genotype.file,head=T,sep=gfile.delimiter,na.strings=c("NA",".","",missing.data))->x
read.table(phenotype.file,head=T,sep=pfile.delimiter,na.strings=c("NA",".","",missing.data))->y
read.table(anno.file,head=T,sep=gfile.delimiter,na.strings=c("NA",".","",missing.data))->m

if (vtype.use[1]!="ALL") m=m[(m[,vtype.col] %in% vtype.use),]
if (covariates=="NONE") covariates=NULL

######################################### collapsing.test()

collapsing.test=function(x,y,m,trait,gene,
ytype,covar=NULL,class_unit,gvid,avid,gsubid,psubid,maf_cutoff=0.01,permun=10000,ped=NULL)

{

### select data for analysis
xs=as.character(m[m[,class_unit]==gene,avid])
if (!is.null(covar)) {
    covars=strsplit(covar,split="[+]")[[1]]
}
else {
    covars=NULL
}

rst=NULL
if (length(xs)>1)
{

if (gsubid=="FIRST_ROW")
{
x[x[,gvid] %in% xs,] -> x
rownames(x)=x[,gvid]
x[,colnames(x)!=gvid]->x
x=as.data.frame(t(x))
x$FIRST_ROW=rownames(x)
y[,psubid] = gsub("-",".",y[,psubid])
}

colnames(x) %in% xs->xsid
colnames(x)[xsid]=paste("V",colnames(x)[xsid],sep="")
xs=colnames(x)[xsid]

y[,colnames(y) %in% c(psubid,trait,covars)]->y
y=y[rowSums(is.na(y))==0,]
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

z=list(design=design,id=id,grp=grp,y=yi,x=xi,ped=ped,ytype=ytype,data=xy[,c(trait,covars,xs)],
trait=trait,variants=xs,covar=covar,covars=covars)
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

date()

try(collapsing.test(x=x,y=y,m=m,
trait=trait,gene=gene,ytype=trait.type,
covar=covariates,
gvid=gfile.vid,
avid=afile.vid,
class_unit=gene.col,
gsubid=gfile.sid,
psubid=pfile.sid,
maf_cutoff=maf.cutoff,
permun=permu.n))->rst

date()

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
WSS=rst$pp["WSS"]
aSum=rst$pp["aSum"]
PWST=rst$pp["PWST"]
SPWST=rst$pp["SPWST"]
SPWST.up=rst$pp["xuplog"]
SPWST.down=rst$pp["xdwlog"]
pp=cbind(Trait,Gene,N,V,MAF,CMC,pCMC,WSS,aSum,PWST,SPWST,SPWST.up,SPWST.down)

write.csv(pp,file=paste(out.dir,"/",trait,"_",gene,".burden.csv",sep=""),quote=F,row.names=F)

p=rst$svt
Variant=rownames(p)
colnames(p)=c("Beta","SE","t","P")
p=signif(p,5)
p=cbind(Variant,p)

write.csv(p,file=paste(out.dir,"/",trait,"_",gene,".single.csv",sep=""),quote=F,row.names=F)

}
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

}


