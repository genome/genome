###########################################################################
# This file contains a collection of R functions for statistical analysis
# of genome-wide data, by Qunyuan Zhang, DSG.  03-17-2008 updated

#Updated 04-08-2010 by NDees for use in the Copy Number Pipeline
###########################################################################
cnn_loh_graph=function(varscan_file,output_file)
{

x<-read.table(file=varscan_file,sep="\t",quote = "");
x$V8 = as.numeric(sub('%','',(x$V8)));
x$V12 = as.numeric(sub('%','',(x$V12)));
pdf(file=output_file,width=10, height=7.5);
plot(x$V2,x$V8,col='blue',type="p",pch=19,xlab=paste("Chr ",x$V1[1]," Coordinates"),ylab="Tumor (Red) & Normal (Blue) Var. Allele Freq",ylim=c(0,100),cex=0.1);
points(x$V2,x$V12,col='red',type="p",pch=19,cex=0.1);
dev.off();

}
###########################################################################
per_sample_cn_graph=function(workdir)
{

mapfile=paste(workdir,'/map.csv',sep='');
CNdir=paste(workdir,'/CN',sep='');
IMGdir=paste(workdir,'/cn_png/',sep='');
dir.create(IMGdir,recursive=T);

read.table(mapfile,header=T)->map;
dir(CNdir)->fs;
for (fi in fs)
  {
    #ffi=paste("/gscmnt/sata181/info/medseq/llin/Ovarian/SNP/hudsonalpha.org_OV.Human1MDuo.1.0.0_analysis/CN", fi, sep="/")
    ffi=paste(CNdir,fi,sep="/")

    log2rdata=ffi
    read.table(log2rdata, header=T)->z
    z=z[,1];

    mapnew=cbind(map,z)
    imgFile=paste(IMGdir,fi,sep="")

    #dl.plotgenome(mapnew, y="z",cutoff=NULL,cutline=0,img=imgFile,yscale=NULL,draw=TRUE,chrom=NULL,mbp=NULL,chr.col="CHR",pos.col="POS",tombp=T)
    plotgenome(mapnew, y="z",cutoff=NULL,cutline=0,img=imgFile,yscale=c(-3,3),draw=TRUE,chrom=NULL,mbp=NULL,chr.col="CHR",pos.col="POS",tombp=T)
  }
}
############################################################################
iscan_per_sample_cn_graph=function(workdir)
{
OUTdir=paste(workdir,'/raw_images/',sep="");
dir.create(OUTdir,recursive=T);

dir(workdir)->fs;
for (fi in fs)
    {
    read.table(fi,header=F)->cn;
    colnames(cn) = c("CHR","POS","cn");
    imgFile=paste(OUTdir,fi,sep="/");

    plotgenome(cn, y="cn",cutoff=NULL,cutline=0,img=imgFile,yscale=c(-3,3),draw=TRUE,chrom=NULL,mbp=NULL,chr.col="CHR",pos.col="POS",tombp=T)
    }
}
########################################################################
run_swt_test=function(workdir)
{

OUTdir=paste(workdir,'/SWT_30_10/',sep="");
dir.create(OUTdir,recursive=T);

mapfile=paste(workdir,"/map.csv",sep="");
read.table(mapfile,header=T)->map;
map = map[,c("CHR","POS")]

dir(paste(workdir,"/CN/",sep=""))->fs
for (fi in fs)
  {
    ffi=paste(workdir,"/CN/", fi, sep="")

    log2file=ffi

    read.table(log2file, header=T)->z
    z = z[,1]

    cbind(map,z)->q
    colnames(q) = c("CHR","POS","cn")

    rFile=paste(OUTdir,fi,".swt",sep="");
    cn.swt(x=q,out.rdata=rFile,wsize=30,wby=10)
  }

}
#######################################################################
iscan_run_swt_test=function(workdir)
{

OUTdir = paste(workdir,'/SWT_30_10/',sep="");
dir.create(OUTdir,recursive=T);
dir(workdir)->fs;
for (fi in fs)
    {
    read.table(fi,header=F)->cn;
    colnames(cn) = c("CHR","POS","cn");
    output_file=paste(OUTdir,fi,".swt",sep="");
    cn.swt(x=cn,out.file=output_file,wsize=30,wby=10);
    }
}
##################################################################
per_sample_swt_graphs=function(workdir)
{

IMGdir=paste(workdir,'/swt_png_30_10/',sep="");
dir.create(IMGdir,recursive=T);

dir(paste(workdir,"/SWT_30_10/", sep=""))->fs
for(fi in fs)
  {
    swtfile = paste(workdir,"SWT_30_10", fi, sep="/")
    load(swtfile)

    z$cn=2*2^z$cn
    graphFile=paste(workdir,"/swt_png_30_10/",fi, sep="")
    plotgenome(z,y="cn", yscale=c(0,6),chr.col="chrom",pos.col="pos1",img=graphFile)

  }
}
##################################################################
iscan_per_sample_swt_graphs=function(workdir)
{
IMGdir=paste(workdir,'/swt_png_30_10/',sep="");
dir.create(IMGdir,recursive=T);
dir(paste(workdir,"/SWT_30_10/", sep=""))->fs;
for(fi in fs)
  {
    swtfile = paste(workdir,"SWT_30_10", fi, sep="/")
    load(swtfile)

    z$cn=2*2^z$cn
    graphFile=paste(workdir,"/swt_png_30_10/",fi, sep="")
    plotgenome(z,y="cn", yscale=c(0,6),chr.col="chrom",pos.col="pos1",img=graphFile)

  }
}
#######################################################
p.correct=function(p)
{
logp=-log10(p)
p.adjust(p=p,method="BH")->fdr
p.adjust(p=p,method="bonferroni")->bonferroni
p=cbind(p,logp,fdr,bonferroni)
p
}
########################################################
g.freq.test=function(x,alt="greater",type="binom")
{
p=numeric(0)
for (i in c(1:nrow(x)))
{
xi=as.numeric(as.character(x[i,]))
if (type=="binom") {pi=binom.test(xi[1],xi[2],xi[3],alternative=alt)$p.value}
if (type=="prop2") {pi=prop.test(x=xi[1:2],n=xi[3:4],alternative=alt)$p.value}
if (type=="fisher") {pi=fisher.test(x=matrix(xi,2,2),alternative=alt)$p.value}
p=c(p,pi)
} 
p.correct(p)->p
invisible(p)
}
############################################################
mutation_freq_test=function(inf,outf=NULL,backgrd=NULL,option=NULL,sep="\t",alt="greater")
{

#data format 1: gene, mutN, coverage 
#data format 2: gene, mutN, coverage, backgroundRate

mut=inf
if (length(inf)==1) {read.table(inf,sep=sep)->mut}

if (!is.null(backgrd))
{

if (is.null(option))
{
expected_n=mut[,3]*backgrd; mut=cbind(mut,backgrd,expected_n)
x=cbind(mut[,2:3],backgrd)
p=g.freq.test(x,alt,type="binom")
}

if (!is.null(option))
{
	if (option=="prop2")
	{
tn=sum(mut[,3])
expected_n=tn*backgrd; mut=cbind(mut,expected_n,tn)
x=mut[,c(2,4,3,5)]
p=g.freq.test(x,alt,type="prop2")
	}
}

}

if (is.null(backgrd))
{

if (!is.null(option))
{
	if (option=="protein")
	{
	freq1=mut[,3]/mut[,5]; freq2=mut[,4]/mut[,6]
	mut=cbind(mut,freq1,freq2); x=mut[,c(3:6)]
	p=g.freq.test(x,alt,type="prop2")
	}

	if (option=="prop2")
	{
expected_n=mut[,3]*mut[,4]; mut=cbind(mut,expected_n)
x=mut[,c(2,5,3,3)]
p=g.freq.test(x,alt,type="prop2")
	}

}

if (is.null(option))
{
expected_n=mut[,3]*mut[,4]; mut=cbind(mut,expected_n)
x=mut[,2:4]
p=g.freq.test(x,alt,type="binom")
}

}

mut=cbind(mut,p)
if (!is.null(outf)) {write.table(mut,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(mut)

}
###################################################################################
mutation_drive_test_bak=function(inf="a data matrix or a file",outf=NULL,alt="greater",sep="")
{
if (length(inf)==1) {read.table(inf,header=T,sep=sep)[,-1]->inf}
inf=inf[,colSums(inf)>0]
inf[inf>1]=1
gene=colnames(inf)
tt=numeric(0)
for (gi in gene)
	{
	inf[,gi]==1->id
	gene==gi->gid
	temp=inf[,gi];x0=sum(temp);n0=length(temp);f0=x0/n0
	temp=inf[id,!gid];x1=sum(temp);n1=nrow(temp)*ncol(temp);f1=x1/n1
	temp=inf[!id,!gid];x2=sum(temp);n2=nrow(temp)*ncol(temp);f2=x2/n2
	temp=cbind(x0,n0,f0,x1,n1,f1,x2,n2,f2)
	tt=rbind(tt,temp)
	}
p=g.freq.test(x=tt[,c("x1","x2","n1","n2")],alt=alt,type="prop2")
tt=cbind(gene,as.data.frame(tt),p)
tt=tt[order(tt[,"p"]),]
if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)
}
###################################################################################
mutation_drive_test=function(inf="a data matrix or a file",outf=NULL,alt="greater",sep="",method="fisher")
{
if (length(inf)==1) {read.table(inf,header=T,sep=sep)[,-1]->inf}
inf[inf>1]=1
gene=colnames(inf)[colSums(inf)>0 & colSums(inf)<nrow(inf)]
tt=numeric(0)
for (gi in gene)
	{
	inf[,gi]==1->id
	gene==gi->gid
	temp=inf[,gi];x0=sum(temp);n0=length(temp);f0=x0/n0
	temp=inf[id,!gid];x1=sum(temp);n1=nrow(temp)*ncol(temp);f1=x1/n1
	temp=inf[!id,!gid];x2=sum(temp);n2=nrow(temp)*ncol(temp);f2=x2/n2
	temp=cbind(x0,n0,f0,x1,n1,f1,x2,n2,f2)
	tt=rbind(tt,temp)
	}
p=g.freq.test(x=tt[,c("x1","x2","n1","n2")],alt=alt,type=method)
tt=cbind(gene,as.data.frame(tt),p)
tt=tt[order(tt[,"p"]),]
if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)
}


###################################################################################

mutation_drive_test.bp=function(inf="a data matrix or a file",bp="gene.bp matrix or file",
outf=NULL,alt="greater",sep="",method="fisher")

{
if (length(inf)==1) {read.table(inf,header=T,sep=sep)->inf; rownames(inf)=inf[,1];inf=inf[,-1]}
if (length(bp)==1) {read.table(bp,header=T,sep=sep)->bp;rownames(bp)=bp[,1];bp=bp[,-1]}
gs=intersect(colnames(inf),colnames(bp))
id=intersect(rownames(inf),rownames(bp))
inf=inf[id,gs];bp=bp[id,gs]

#inf[inf>1]=1

gene=gs[colSums(inf)>0]
tt=numeric(0)
for (gi in gene)
	{
	inf[,gi]>0->id  # sample id
	gene==gi->gid
	temp=inf[,gi];x0=sum(temp);temp=bp[,gi];n0=sum(temp);f0=x0/n0
	temp=inf[id,!gid];x1=sum(temp);temp=bp[id,!gid];n1=sum(temp);f1=x1/n1
	temp=inf[!id,!gid];x2=sum(temp);temp=bp[!id,!gid];n2=sum(temp);f2=x2/n2
	temp=cbind(x0,n0,f0,x1,n1,f1,x2,n2,f2)
	tt=rbind(tt,temp)
	}
p=g.freq.test(x=tt[,c("x1","x2","n1","n2")],alt=alt,type=method)
tt=cbind(gene,as.data.frame(tt),p)
tt=tt[order(tt[,"p"]),]
if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)
}

###########################################
pathway_test=function
(x,gcolname="gene",lev=c("lev1","lev2","lev3"),yn=c("yn1","yn2","yn3"),cutoff=0.5,outf=NULL,sep="")
{
if (length(x)==1) {read.table(x,header=T,sep=sep)->tt} else {tt=x;rm(x)}
lev=intersect(lev,colnames(tt))
tttb=numeric(0)
for (yni in yn)
{
unique(tt[,c(gcolname,yni)])[,yni]->temp
up_all=sum(temp>=cutoff,na.rm=T);down_all=sum(temp<cutoff,na.rm=T);all=up_all+down_all
for (levi in lev)
{ 
pathnames=levels(as.factor(tt[,levi]))
ttb=numeric(0)
for (pathi in pathnames)
{
tt[(tt[,levi]==pathi),c(gcolname,yni)]->a
all_in=length(unique(a[,gcolname]));all_out=all-all_in
unique(a)[,yni]->temp
up_in=sum(temp>=cutoff,na.rm=T);down_in=sum(temp<cutoff,na.rm=T)
up_out=up_all-up_in;down_out=down_all-down_in;
level=levi;pathway=pathi;index=yni
tb=cbind(index,level,pathway,all_in,all_out,up_in,down_in,up_out,down_out)
ttb=rbind(ttb,tb)
}
g.freq.test(x=ttb[,c("up_in","down_in","up_out","down_out")],alt="two.sided",type="fisher")->p
ttb=as.data.frame(ttb);ttb=cbind(ttb,p)
ttb=ttb[order(ttb[,"p"]),]
tttb=rbind(tttb,ttb)
}
}
if (!is.null(outf)) {write.table(tttb,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tttb)
}

############################################################
mutation_pathway_test=function
(inf_mut=NULL,inf_path=NULL,outf=NULL,sigmut=1,backgrd=NULL,alpha=0.05,
gcolname="gene",lev=3,mutinfo=c("mutn","bp"))
{

if (length(inf_mut)>1) {tt=inf_mut;rm(inf_mut)} else merge_mut_pathway_file(inf_mut,inf_path)->tt 
if (is.null(backgrd)) {tt[,mutinfo[1]]>=sigmut->ismut} else
{
ismut=cbind(tt[,mutinfo],backgrd)
ismut=g.freq.test(ismut,alt="greater",type="binom")
ismut[,1]<alpha -> ismut
}
ismut=as.numeric(ismut)
tt=cbind(tt,ismut)

if (class(lev)=="numeric") {lev=paste("level",c(1:lev),sep="_")} 
lev=intersect(lev,colnames(tt))

pathway_test(tt[,c(gcolname,lev,"ismut")],gcolname,lev,"ismut")->gpt
gpt=gpt[,-1]; colnames(gpt)[9:12]=paste("g",colnames(gpt)[9:12],sep=".")

#-----------------------------------------------------------------------

colSums(unique(tt[,c(gcolname,mutinfo)])[,mutinfo])->mutsum
mutn0=mutsum[1];bp0=mutsum[2];gene_all=length(unique(tt[,gcolname]))

tttb=numeric(0)
for (level in lev)
{ 
ttb=numeric(0)
pathnames=levels(as.factor(tt[,level]))
for (pathi in pathnames)
{
a=tt[tt[,level]==pathi,c(gcolname,mutinfo)]
gene_in=length(unique(a[,gcolname]))
gene_out=gene_all-gene_in
colSums(unique(a)[,mutinfo])->pisum
mut_in=pisum[1];bp_in=pisum[2];freq_in=mut_in/bp_in
mut_out=mutn0-mut_in;bp_out=bp0-bp_in;freq_out=mut_out/bp_out
pathway=pathi
tb=cbind(level,pathway,gene_in,gene_out,mut_in,bp_in,mut_out,bp_out,freq_in,freq_out)
ttb=rbind(ttb,tb)
}

g.freq.test(x=ttb[,c("mut_in","mut_out","bp_in","bp_out")],alt="two.sided",type="prop2")->p
ttb=cbind(ttb,p)
if (!is.null(backgrd))
{
g.freq.test(x=cbind(ttb[,c("mut_in","bp_in")],backgrd),alt="greater",type="binom")->p
colnames(p)=paste("b",colnames(p),sep=".")
ttb=cbind(ttb,p)
}
ttb=as.data.frame(ttb)

tttb=rbind(tttb,ttb)
}

tttb=merge(tttb,gpt,all=T)
tttb=tttb[order(tttb[,"p"]),]
tttb[,"pathway"]=gsub(","," ",as.character(tttb[,"pathway"]))
tttb[,-(1:2)]=as.numeric.frame(tttb[,-(1:2)])
if (!is.null(outf)) {write.table(tttb,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tttb)

}
###############################

merge_mut_pathway_file=function(mutf,pathf)
{
read.table(pathf,sep="\t",stringsAsFactors=F)->pa
colnames(pa)=c("id","gene",paste("level",c(1:(ncol(pa)-2)),sep="_"))
read.table(mutf,sep="\t",stringsAsFactors=F)->mut
colnames(mut)=c("gene","mutn","bp")
x=merge(mut,pa,by.x="gene",by.y="gene")
x
}

#################################

cor2test_bak =function(y,x=NULL,method="cor",cutoff=1,sep="",outf=NULL)
{

if (length(x)==1) {read.table(x,header=T,sep=sep)->x}
if (length(y)==1) {read.table(y,header=T,sep=sep)->y}
if (is.null(x)) {x=y}
colnames(y)[1]="id";colnames(x)[1]="id"

tt=character(0)
for (vi in colnames(x)[-1])
{
for (vj in colnames(y)[-1])
{
tx=x[,c("id",vi)];tx=tx[!is.na(tx[,vi]),];tx=tx[!duplicated(tx[,"id"]),]
ty=y[,c("id",vj)];ty=ty[!is.na(ty[,vj]),];ty=ty[!duplicated(ty[,"id"]),]
xy=merge(tx,ty,by.x="id",by.y="id")
n=length(xy[,"id"])

if (n>5)
{
tx=xy[,2];ty=xy[,3];p=NA;s=NA

if(method=="cor")
{
tst=cor.test(tx,ty)
s=tst$est
p=tst$p.value
}

if(method=="chisq")
{
levx=length(levels(factor(tx)))
levy=length(levels(factor(ty)))
if(levx>1 & levx<15 & levy>1 & levy<15)
{
tst=chisq.test(tx,ty)
s=tst$stat
p=tst$p.value
}
}

if(method=="fisher")
{
levx=length(levels(factor(tx)))
levy=length(levels(factor(ty)))
if(levx>1 & levx<5 & levy>1 & levy<5)
{
tst=fisher.test(tx,ty)
s=tst$p.value
p=tst$p.value
}
}

if(method=="anova")
{
levx=length(levels(factor(tx)))
#levy=length(levels(factor(ty)))
if(levx>1 & levx<10)
{
tst=summary(aov(ty~tx,as.data.frame(cbind(tx,ty))))
s=tst[[1]]$F[1]
p=tst[[1]]$Pr[1]
}
}

t=c(vi,vj,method,n,s,p)
tt=rbind(tt,t)

} #end if n>3
} #end vj
} #end vi

rownames(tt)=NULL
colnames(tt)=c("x","y","method","n","s","p")
tt=as.data.frame(tt)
tt[,"s"]=as.character(tt[,"s"]);tt[,"s"]=as.numeric(tt[,"s"])
tt[,"p"]=as.character(tt[,"p"]);tt[,"p"]=as.numeric(tt[,"p"])
fdr=p.adjust(tt[,"p"],method="fdr")
bon=p.adjust(tt[,"p"],method="bon")
tt=cbind(tt,fdr,bon)
tt=tt[tt[,"p"]<=cutoff,]
tt=tt[order(tt[,"p"]),]

if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)

} #end cortest_2var

################################

cor2test =function(y,x=NULL,method="cor",cutoff=1,sep="",outf=NULL)
{


if (!is.null(x))
{

if (length(x)==1) {read.table(x,header=T,sep=sep)->x}
if (length(y)==1) {read.table(y,header=T,sep=sep)->y}
colnames(y)[1]="id";colnames(x)[1]="id"
tt=character(0)
for (vi in colnames(x)[-1])
{
for (vj in colnames(y)[-1])
{
tx=x[,c("id",vi)];tx=tx[!is.na(tx[,vi]),];tx=tx[!duplicated(tx[,"id"]),]
ty=y[,c("id",vj)];ty=ty[!is.na(ty[,vj]),];ty=ty[!duplicated(ty[,"id"]),]
xy=merge(tx,ty,by.x="id",by.y="id")
tx=xy[,2];ty=xy[,3]
n=length(xy[,"id"])
rst=try(cor2(ty,tx,method))
if (class(rst)=="try-error") {p=NA} else {p=rst[1];s=rst[2]}
if (!is.na(p)) { t=c(vi,vj,method,n,s,p); tt=rbind(tt,t) }
} #end vj
} #end vi

rownames(tt)=NULL
colnames(tt)=c("x","y","method","n","s","p")
tt=as.data.frame(tt)
tt[,"s"]=as.character(tt[,"s"]);tt[,"s"]=as.numeric(tt[,"s"])
tt[,"p"]=as.character(tt[,"p"]);tt[,"p"]=as.numeric(tt[,"p"])
fdr=p.adjust(tt[,"p"],method="fdr")
bon=p.adjust(tt[,"p"],method="bon")
tt=cbind(tt,fdr,bon)
tt=tt[tt[,"p"]<=cutoff,]
tt=tt[order(tt[,"p"]),]
}



if (is.null(x))
{

if (length(y)==1) {read.table(y,header=T,sep=sep)->y}
x=y; nxy=ncol(y)-1
colnames(y)[1]="id";colnames(x)[1]="id"
tt=character(0)
for (i in c(1:(nxy-1)))
{
for (j in c((i+1):nxy))
{

vi=colnames(x)[-1][i]
vj=colnames(y)[-1][j]

tx=x[,c("id",vi)];tx=tx[!is.na(tx[,vi]),];tx=tx[!duplicated(tx[,"id"]),]
ty=y[,c("id",vj)];ty=ty[!is.na(ty[,vj]),];ty=ty[!duplicated(ty[,"id"]),]
xy=merge(tx,ty,by.x="id",by.y="id")
tx=xy[,2];ty=xy[,3]
n=length(xy[,"id"])
rst=try(cor2(ty,tx,method))
if (class(rst)=="try-error") {p=NA} else {p=rst[1];s=rst[2]}
if (!is.na(p)) { t=c(vi,vj,method,n,s,p); tt=rbind(tt,t) }
} #end vj
} #end vi

rownames(tt)=NULL
colnames(tt)=c("x","y","method","n","s","p")
tt=as.data.frame(tt)
tt[,"s"]=as.character(tt[,"s"]);tt[,"s"]=as.numeric(tt[,"s"])
tt[,"p"]=as.character(tt[,"p"]);tt[,"p"]=as.numeric(tt[,"p"])
fdr=p.adjust(tt[,"p"],method="fdr")
bon=p.adjust(tt[,"p"],method="bon")
tt=cbind(tt,fdr,bon)
tt=tt[tt[,"p"]<=cutoff,]
tt=tt[order(tt[,"p"]),]
}


if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)

} #end cor2test

################################

cor2test.bak =function(y,x=NULL,method="cor",cutoff=1,sep="",outf=NULL)
{

if (length(x)==1) {read.table(x,header=T,sep=sep)->x}
if (length(y)==1) {read.table(y,header=T,sep=sep)->y}
if (is.null(x)) {x=y}
colnames(y)[1]="id";colnames(x)[1]="id"

tt=character(0)
for (vi in colnames(x)[-1])
{
for (vj in colnames(y)[-1])
{
tx=x[,c("id",vi)];tx=tx[!is.na(tx[,vi]),];tx=tx[!duplicated(tx[,"id"]),]
ty=y[,c("id",vj)];ty=ty[!is.na(ty[,vj]),];ty=ty[!duplicated(ty[,"id"]),]
xy=merge(tx,ty,by.x="id",by.y="id")
tx=xy[,2];ty=xy[,3]
n=length(xy[,"id"])
rst=try(cor2(ty,tx,method))
if (class(rst)=="try-error") {p=NA} else {p=rst[1];s=rst[2]}
if (!is.na(p)) { t=c(vi,vj,method,n,s,p); tt=rbind(tt,t) }
} #end vj
} #end vi

rownames(tt)=NULL
colnames(tt)=c("x","y","method","n","s","p")
tt=as.data.frame(tt)
tt[,"s"]=as.character(tt[,"s"]);tt[,"s"]=as.numeric(tt[,"s"])
tt[,"p"]=as.character(tt[,"p"]);tt[,"p"]=as.numeric(tt[,"p"])
fdr=p.adjust(tt[,"p"],method="fdr")
bon=p.adjust(tt[,"p"],method="bon")
tt=cbind(tt,fdr,bon)
tt=tt[tt[,"p"]<=cutoff,]
tt=tt[order(tt[,"p"]),]

if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)

} #end cor2test

#################################

cor2test0 =function(y,x=NULL,method="cor",cutoff=1,sep="",outf=NULL)
{

if (length(x)==1) {read.table(x,header=T,sep=sep)->x}
if (length(y)==1) {read.table(y,header=T,sep=sep)->y}
if (is.null(x)) {x=y}
colnames(y)[1]="id";colnames(x)[1]="id"

tt=character(0)
for (vi in colnames(x)[-1])
{
for (vj in colnames(y)[-1])
{
tx=x[,c("id",vi)];tx=tx[!is.na(tx[,vi]),];tx=tx[!duplicated(tx[,"id"]),]
ty=y[,c("id",vj)];ty=ty[!is.na(ty[,vj]),];ty=ty[!duplicated(ty[,"id"]),]
xy=merge(tx,ty,by.x="id",by.y="id")
tx=xy[,2];ty=xy[,3]
n=length(xy[,"id"])
rst=try(cor2(ty,tx,method))
if (class(rst)=="try-error") {p=NA} else {p=rst[1];s=rst[2]}
if (!is.na(p)) { t=c(vi,vj,method,n,s,p); tt=rbind(tt,t) }
} #end vj
} #end vi

rownames(tt)=NULL
colnames(tt)=c("x","y","method","n","s","p")
tt=as.data.frame(tt)
tt[,"s"]=as.character(tt[,"s"]);tt[,"s"]=as.numeric(tt[,"s"])
tt[,"p"]=as.character(tt[,"p"]);tt[,"p"]=as.numeric(tt[,"p"])
fdr=p.adjust(tt[,"p"],method="fdr")
bon=p.adjust(tt[,"p"],method="bon")
tt=cbind(tt,fdr,bon)
tt=tt[tt[,"p"]<=cutoff,]
tt=tt[order(tt[,"p"]),]

if (!is.null(outf)) {write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",")}
invisible(tt)

} #end cor2test0

################################

cor2=function(ty,tx,method)
{

id=intersect(!is.na(ty),!is.na(tx))
ty=ty[id];tx=tx[id]

if(method=="cor")
{
tst=cor.test(tx,ty)
s=tst$est
p=tst$p.value
}

if(method=="chisq")
{
tst=chisq.test(tx,ty)
s=tst$stat
p=tst$p.value
}

if(method=="fisher")
{
tst=fisher.test(tx,ty)
s=tst$p.value
p=tst$p.value
}

if(method=="anova")
{
tst=summary(aov(ty~tx,as.data.frame(cbind(tx,ty))))
s=tst[[1]]$F[1]
p=tst[[1]]$Pr[1]
}

tt=c(p,s)
tt

}

######################################################

cn.log2ratio=function(sampleinfo,id=NULL,datadir=NULL,outdir="_log2ratio",sep="",des="tumor/normal pair")
{
if (is.null(datadir)) datadir=getwd()
x=sampleinfo; if (length(x)==1) read.table(paste(datadir,x,sep="/"),header=T,sep=sep,stringsAsFactors=F)->x

x[,2]=paste(datadir,x[,2],sep="/")
x[,3]=paste(datadir,x[,3],sep="/")

if (!file.exists(outdir)) dir.create(outdir)
ids=unique(x[,1])
if (is.null(id)) id=c(1:length(ids))
for (i in id)
{
if (i<=length(ids))
{
xi=x[x[,1]==ids[i],];n=nrow(xi);tum=0;nom=0
for (j in c(1:n))
{
as.numeric(read.table(xi[j,2],header=T)[,1])->tumj;tum=tum+tumj
as.numeric(read.table(xi[j,3],header=T)[,1])->nomj;nom=nom+nomj
}
z=log2(tum/n)-log2(nom/n)
z=as.numeric(format(z,digits=4))
datatype="log2ratio"
rawdata=xi
save(z,datatype,rawdata,des,file=paste(outdir,"/log2ratio","_",ids[i],sep=""),compress=T)
}
}  
}


#######################################################
cn.ave = function(datadir,max.cn=10,sd.range=5,cut.a=2,cut.d=-2)
{
fs=dir(datadir)
fs=paste(datadir,fs,sep="/")
tt=0
for(fi in fs)
    {
load(fi)
m0=mean(z,na.rm=T)
sd0=sd(z,na.rm=T)
m=mean(z[abs((z-m0)/sd0)<sd.range],na.rm=T)
sd=sd(z[abs((z-m0)/sd0)<sd.range],na.rm=T)
sdi=(z-m)/sd
a=integer(length(z));d=a
a[sdi>cut.a]=1;d[sdi<cut.d]=1
z = 2 * 2^z; z[z>max.cn]=max.cn; z[is.na(z)]=m 
z=cbind(z,a,d)	
tt=tt + z
    }
    tt = tt/length(fs);colnames(tt)[1]="cn"
	invisible(tt)
}

########################################################
cn.ave_map = function(datadir,map,max.cn=10,sd.range=5,cut.a=2,cut.d=-2,out.file=NULL,out.rdata=NULL)
{
cn.ave(datadir,max.cn,sd.range,cut.a,cut.d)->cn; cn=as.numeric.frame(cn,fix=4)
read.table(map,header=T)->z
cbind(z,cn)->z
if(!is.null(out.file)) write.table(z,file=out.file,quote=F,row.names=F,sep="\t")
if (!is.null(out.rdata)) save(z,map,file=out.rdata,compress=T)
invisible(z)
}


######################################################

cn.swt=function(x,wsize=30,wby=1,sd.range=5,chr="all",chr.out=FALSE,
chr.col="CHR",pos.col="POS",cn.col="cn",
out.file=NULL,out.rdata=NULL,sep="",draw=F)
{
if (length(x)==1) read.table(x,header=T,sep=sep,stringsAsFactor=F)->x
x=x[,c(chr.col,pos.col,cn.col)]
x[,chr.col]=as.factor(x[,chr.col])
if (chr[1]=="all") chr=levels(x[,chr.col]) else chr=intersect(levels(x[,chr.col]),chr)
x=x[x[,chr.col] %in% chr,]
z=numeric(0)

for (chri in chr)
{
xchri=x[x[,chr.col]==chri,]; nchri=nrow(xchri)
if (nchri>=wsize)
{ 
xchri=xchri[order(xchri[,pos.col]),]
for (wi in seq(1,nchri-wsize+1,by=wby))
{
w = xchri[wi:(wi+wsize-1),]
chrom=chri;pos1=w[1,pos.col];pos2=w[nrow(w),pos.col];cn=mean(w[,cn.col],na.rm=T)
temp=c(chrom,pos1,pos2,cn)
z=rbind(z,temp)
}

if (chr.out==TRUE)
{
colnames(z)=c("chrom","pos1","pos2",cn.col)
write.table(as.numeric.frame(z[z[,1]==chri,],cn.col,fix=4),
file=paste("swt_chrom",chri,".csv",sep=""), quote=F, row.names=F,sep=",")
}

}
}

colnames(z)=c("chrom","pos1","pos2",cn.col)
rownames(z)=NULL
z=as.numeric.frame(z,cn.col,fix=4)
z=as.numeric.frame(z,c("pos1","pos2"),fix=10)
for (cn.coli in cn.col)
{
gvalues.test(z[,cn.coli],sd.range)->p;
colnames(p)=paste(cn.coli,colnames(p),sep=".")
p=as.numeric.frame(p,fix=3)
z=cbind(z,p)
}
if (!is.null(out.file)) write.table(z,file=out.file,quote=FALSE,row.names=FALSE,sep="\t")
if (!is.null(out.rdata)) save(wsize,wby,chr,z,file=out.rdata,compress=T)
invisible(z)
}

###############################################

gvalues.test=function(x,sd.range=5)
{
x=as.numeric(as.character(x))
sd=(x-mean(x[abs(sd(x))<sd.range]))/sqrt(var(x[abs(sd(x))<sd.range]))
p=2*pnorm(-abs(sd));fdr=p.adjust(p,method="BH")
tt=cbind(sd,p,fdr)
invisible(tt)
}

###############################################
as.numeric.frame=function(x,col=NULL,fix=5)
{
x=as.data.frame(x)
if (is.null(col)) col=colnames(x)
for (coli in col) x[,coli]=as.numeric(format(as.numeric(as.character(x[,coli])),digits=fix))
x
}
##############################################

cn.ave_gene = function(gene.snp.file,snp.map.file,log2ratio.dir,max.cn,out.file=NULL,out.rdata=NULL)
{
x=gene.snp.file
map=snp.map.file

if (!file.exists(output.dir)) dir.create(output.dir)
setwd(output.dir)

if (class(x)=="character")
{read.table(x,header=T,stringsAsFactors=F)[,c(probe.col,gene.col),]->x}
x=x[x$HUGO!="--",]

if (class(map)=="character")
{read.table(map,header=T,stringsAsFactors=F)[,c(probe.col)]->map}
map=as.data.frame(map);colnames(map)=probe.col

fs=paste(log2ratio.dir,dir(log2ratio.dir),sep="/")

tt=NULL;ID=NULL
for(fi in fs)
    {
load(fi); z=2*2^z;z[z>max.cn]=max.cn
z=cbind(map,z)
merge(x,z,by.x=probe.col,by.y=probe.col)->x
tapply (x[,"z"],x[,gene.col],mean,na.rm=T)->temp
x=x[,c(probe.col,gene.col)]
tt=cbind(tt,temp)
strsplit(fi,"_")[[1]]->idi; idi=idi[length(idi)]
ID=c(ID,idi)
    }

colnames(tt)=ID; tt=as.numeric.frame(tt,fix=3)   
if(!is.null(out.file)) write.table(tt,file=out.file,quote=F,sep="\t")
if (!is.null(out.rdata)) save(tt,file=out.rdata,compress=T)
invisible(tt)
}

################################################################

cn.gene.test.all=function(x,out.file=NULL,out.rdata=NULL,cn.lim=4,sep="")
{
if (length(x)==1) read.table(x,header=T,sep=sep,stringsAsFactor=F)->x
x=t(x)
m0=mean(x[x<cn.lim],na.rm=T)
m=colMeans(x,na.rm=T)
s=sd(x,na.rm=T)
n=colSums(!is.na(x))
t=(m-m0)/s*sqrt(n)
p=pt(-abs(t),n-1)*2
p=p.correct(p)
z=cbind(m0,m,s,n,t,p)
z=as.numeric.frame(z)
if (!is.null(out.file)) write.csv(z,file=out.file,quote=FALSE)
des="gene CN test over entire samples"
if (!is.null(out.rdata)) save(z,des,file=out.rdata,compress=T)
invisible(z)
}



################################################################

cn.smooth.indv=function(datadir=NULL,id=NULL,outdir="_smooth.indv")
{
if (is.null(datadir)) datadir=getwd()
fs0=dir(datadir)
fs=paste(datadir,fs0,sep="/")
if (!file.exists(outdir)) dir.create(outdir)
if (is.null(id)) id=c(1:length(fs))
for (i in id)
{
if (i>0 & i<=length(fs))
{
load(fs[i])
z=z[1:5]
datatype="CN"
rawdata=fs[i]
save(z,datatype,rawdata,file=paste(outdir,"/",fs0[i],".glad",sep=""),compress=T)
}
}  

}

#####################################################################################
plotgenome = function (tt, y="p",cutoff=NULL,cutline=2,img=NULL,yscale=NULL,draw=TRUE,ltype="p",
chrom=NULL,mbp=NULL,chr.col="chromosome",pos.col="position",tombp=T)
{
colnames(tt)[grep(chr.col,colnames(tt))]="chr"
colnames(tt)[grep(pos.col,colnames(tt))]="mbp"
if (tombp) tt$mbp=tt$mbp/1000000
if (is.null(chrom)) {chrom=levels(as.factor(tt[,"chr"]));chrom=intersect(c(1:24,"X","Y","x","y"),chrom)}
if (length(mbp)==2) {tt=tt[(tt[,"mbp"]>=mbp[1] & tt[,"mbp"]<=mbp[2]),]}
if (!is.null(cutoff)) {tt[,y][tt[,y]<cutoff]=NA}
if (draw==TRUE)
{
#if (length(img)>0) {png(paste(img,"png",sep="."),1200,800)}
if (length(img)>0) {bitmap(file=paste(img,"png",sep="."),height=11,width=16, res=200)}

length(chrom)->chrn
if (chrn==1){mr=1;mc=1}
if (chrn==2){mr=2;mc=1}
if (chrn==3){mr=2;mc=2}
if (chrn==4){mr=2;mc=2}
if (chrn==5){mr=2;mc=3}
if (chrn==6){mr=2;mc=3}
if (chrn==7){mr=3;mc=3}
if (chrn==8){mr=3;mc=3}
if (chrn==9){mr=3;mc=3}
if (chrn==10){mr=3;mc=4}
if (chrn==11){mr=3;mc=4}
if (chrn==12){mr=3;mc=4}
if (chrn==13){mr=3;mc=5}
if (chrn==14){mr=3;mc=5}
if (chrn==15){mr=3;mc=5}
if (chrn==16){mr=4;mc=4}
if (chrn==17){mr=4;mc=5}
if (chrn==18){mr=4;mc=5}
if (chrn==19){mr=4;mc=5}
if (chrn==20){mr=4;mc=5}
if (chrn==21){mr=4;mc=6}
if (chrn==22){mr=4;mc=6}
if (chrn==23){mr=4;mc=6}
if (chrn==24){mr=4;mc=6}
if (chrn>24){mr=5;mc=6}

par(mfrow=c(mr,mc))
for (chr in chrom )
{
tt[,"chr"]==chr->ch
tl=paste("Chrom.",chr,sep=" ")
matplot(tt[ch,"mbp"],tt[ch,y],pch=".",main=tl,xlab="Mbp",ylab=y,ylim=yscale,type=ltype)
for (ct in cutline) {abline(ct,0,col="red")}
}#chr
if (length(img)>0) {dev.off()}
}#draw

invisible(tt)
}
###############################################################################
seq.vjoint=function(x)
{
tt=NULL
for (i in 1:length(x)) {tt=c(tt,rep(names(x[i]),x[i]))}
tt
}

###############################################################################
seq.mjoint=function(x)
{
tt=NULL
for (i in 1:nrow(x))
{
seq.vjoint(x[i,])->xi
if (length(xi)>0) {tt=rbind(tt,cbind(i,xi))}
}
tt 
}
###############################################################################
prob.concur=function(x)
{
x>0->x
cols=colnames(x)
#
read.table("/gscuser/lding/TSP_manuscript/Mutation_Matrix/TSP_Lung_Adenocarcinoma_Mutation_Table.071408.gene_6mut.matrix.csv",header=T)->subx
cols=intersect(cols,colnames(subx))
#
cols=cols[order(cols)]; n=length(cols)

tt=NULL
for (i in 1:(n-1))
{
for (j in (i+1):n)
{
vi=cols[i];vj=cols[j]
sum(x[,vi] & x[,vj])->nand
sum(x[,vi] != x[,vj])->nexc
ni=sum(x[,vi]);nj=sum(x[,vj])
temp=cbind(vi,vj,ni,nj,nand,nexc)
tt=rbind(tt,temp)
}
}
tt=as.data.frame(tt,stringsAsFactors=F)
tt$nand=as.numeric(tt[,"nand"])
tt$nexc=as.numeric(tt[,"nexc"])
tt
}
###############################################################################
# FUNCTION without counting gene-specific mut freq., sample by sample
matrix.sample=function(x,keep.freq=T)
{
colSums(x)->pool
rowSums(x)->freq;freq=freq[order(freq,decreasing=T)]
tt=NULL
for (i in c(1:length(freq)))
{
if (length(pool)>1) sample(names(pool),freq[i])->xi 
if (length(pool)==1) names(pool)->xi
pool[xi]=pool[xi]-1;pool=pool[pool>0]
tt=rbind(tt,cbind(i,xi))
}
tt=table(tt[,1],tt[,2])
tt
}
###############################################################################
# FUNCTION counting gene-specific mut freq., sample by sample
matrix.sample.v1=function(x,keep.freq=T)
{
colSums(x)->pool
rowSums(x)->freq;freq=freq[order(freq,decreasing=T)]
tt=NULL
for (i in c(1:length(freq)))
{
pooli=seq.vjoint(pool)
xi=NULL
for (j in 1:freq[i])
{
if (length(pooli)>1) sample(pooli,1)->temp 
if (length(pooli)==1) pooli->temp
xi=c(xi,temp)
pooli=pooli[pooli!=temp]
}
tt=rbind(tt,cbind(i,xi))
pool[xi]=pool[xi]-1;pool=pool[pool>0]
}
tt=table(tt[,1],tt[,2])
tt
}
###############################################################################
# FUNCTION counting gene-specific mut freq., random sample 
seq.matrix.sample.v2=function(x,keep.freq=T)
{
x=as.data.frame(x);x[,2]=as.character(x[,2]);x[,1]=as.character(x[,1])
x1=x;x1[,2]=sample(x[,2]);x1$uni=NA
x2=unique(x1);x1[rownames(x2),3]=as.character(x2[,2])
x1[is.na(x1$uni),]->x2
genes=sample(x2[,2])

for (gi in genes[1:10])
{
for (si in rownames(x2))
{
if (!(gi %in% x1[x1[,"i"]==x2[si,"i"],"uni"])) 
{x1[si,"uni"]=gi; x2=x2[rownames(x2)!=si,]; break}
}
}
x2
x1[is.na(x1$uni),]->x2
x2

tt=table(x[,1],x[,2])
tt
}
###############################################################################
mut_cor_permu_test=function(x=NLL,n.permu=100,seed=NULL,
mut.file="/gscuser/qzhang/gstat/example_data/WU_broad_4mutations.matrix.csv",
out.file=NULL,out.rdata=NULL)

{
if (!is.null(mut.file)) read.table(mut.file,header=T)->x
x=x[,-1];x[x>0]=1
x=x[rowSums(x)>0,colSums(x)>0]
prob.concur(x)->pc0
#seq.mjoint(x)->seq.x

if (!is.null(seed)) set.seed(seed)

pp=0;en=0
for (i in 1:n.permu)
{
#xi=cbind(sample(x0[,1]),sample(x0[,2]))

#xi=matrix.sample.v1(x,keep.freq=T)

doit=1; while(doit==1) {xi=matrix.sample.v1(x); if (max(xi)==1) doit=0} 

#xi=seq.matrix.sample.v2(seq.x,keep.freq=T)

prob.concur(xi)->pci
#head(pci)
#head(merge(pc0,pci,by.x=c("vi","vj"),by.y=c("vi","vj")))
en=en+pci[,c("nand","nexc")]
as.numeric(pci$nand>=pc0$nand)->pand
as.numeric(pci$nexc>=pc0$nexc)->pexc
pp=pp+cbind(pand,pexc)
}
en=en/n.permu;pp=pp/n.permu
pci[,c("nand","nexc")]=en
pci=cbind(pci,pp)
tt=merge(pc0,pci,by.x=c("vi","vj"),by.y=c("vi","vj"))

if (!is.null(out.rdata)) save(tt,file=out.rdata,compress=T)
if (!is.null(out.file)) write.csv(tt,file=out.file,quote=F)

invisible (tt)

}
####################################################################################
