##############################################################
# CN SIMULATION FOR SMDS ...
##############################################################

Simulate_Regions=function(n,regions,cor=FALSE,viberate=NULL,vary=T)
{
# background intensity
row_n=n
col_n=regions["background","end"]
mu=regions["background","mean"]
sigma=regions["background","sd"]
simu_data=matrix(0,row_n,col_n)
simu_data[,]=rnorm(row_n*col_n,mean=mu,sd=sigma)

# non-recurrent regions
regions["nonrec","start"]->lamda
if (lamda>0)
{
nonrec_n=as.integer(regions["nonrec","freq"]*row_n)
nonrec_len=rpois(nonrec_n,lamda)+1
start=as.integer(runif(nonrec_n)*(col_n-nonrec_len))
start[start==0]=1
sampid=sample(1:row_n,nonrec_n)
sigma=regions["nonrec","sd"]
mu0=regions["nonrec","mean"]
for(i in c(1:length(sampid)))
{
if (is.null(viberate)) mu=mu0 else mu=cn.viberate(viberate) 
simu_data[sampid[i],start[i]:(start[i]+nonrec_len[i]-1)]=rnorm(nonrec_len[i],mean=mu,sd=sigma)
}

}

# recurrent regions
for (region_i in rownames(regions)[-(1:2)])
{
start=regions[region_i,"start"]
end=regions[region_i,"end"]
common_n=as.integer(regions[region_i,"freq"]*row_n)
mu0=regions[region_i,"mean"]
sigma=regions[region_i,"sd"]
for(i in c(1:common_n))
{
if (is.null(viberate)) mu=mu0 else mu=cn.viberate(viberate)
if (vary==F) {endi=end;starti=start}
if (vary==T)
{mid=(end-start)/2+start;vi=rpois(1,end-start)
endi=mid+vi/2;starti=mid-vi/2}
 
simu_region=rnorm((endi-starti+1),mean=mu,sd=sigma)
if (cor==TRUE) {simu_data[i,starti:endi]=simu_region}
if (cor==FALSE) {simu_data[sample(1:row_n,1),starti:endi]=simu_region}
}
}
simu_data
}


cn.viberate=function(f)
{
f=f/sum(f)
m=log2((c(1:length(f))-1)/2);m[1]=-3
f2=f;for(i in c(1:length(f))) f2[i]=sum(f[1:i])
f1=c(0,f2[-length(f2)]) 
x=runif(1); mx=m[x>f1 & x<f2]
mx 
}

##############################################################
# CMDS ...
##############################################################

##############################################################
cmds=function(x,wsize=10,wstep=1,ztrans=T,alpha=0.05,m.cut=0)
{
Diagonalize_Regions(x,wsize,wstep,ztrans)[,"z"]->temp
temp=cbind(c(1:length(temp)),temp)
Smooth_Regions(temp,norm=T,draw=F)-> smo
smo$adj_region=smo$region
Test_Regions(smo,alpha,m.cut)->smo
Adjust_Sig_Rigions(smo,tail=wsize,ignore=wsize)->smo
invisible(smo)
}

#########################################################

Smooth_Regions=function(c,norm=T,draw=F)
{
x0=c[,2]
if (norm==TRUE){c[,2]=(c[,2]-mean(c[,2]))/sd(c[,2])}
c=as.data.frame(c)
cgh=cbind(c(1:dim(c)[1]),c,1)
colnames(cgh)=c("PosOrder","PosBase","LogRatio","Chromosome")
cgh<- as.profileCGH(cgh)

 res <- glad(cgh, mediancenter = FALSE, smoothfunc = "lawsglad",
 bandwidth = 10, round = 1.5, model = "Gaussian", lkern = "Exponential",
 qlambda = 0.99, base = FALSE, lambdabreak = 8, lambdacluster = 8,
 lambdaclusterGen = 40, type = "tricubic", param = c(d = 6),
 alpha = 0.001, msize = 5, method = "centroid", nmax = 8,
 verbose = FALSE)

res=res$pro[,c("PosOrder","PosBase","LogRatio","Smoothing","Region")]
colnames(res)=c("order","pos","x","smth","region")
res=cbind(res,x0)

#plotProfile(res, unit = 3, Bkp = TRUE, labels = FALSE, Smoothing = "Smoothing", plotband = FALSE)
if (draw==TRUE)
{
plot(res[,"x"],type="l")
lines(res[,"smth"],col="Blue",lwd=2)
abline(0,0,col="Red",lty=3)
abline(2,0,col="Red",lty=3)
abline(-2,0,col="Red",lty=3)
}
res
}
 
##############################################################

Test_Regions=function(x,alpha=0.05,m.cut=0)
{
tapply(x$x,x$adj_region,mean)->m
tapply(x$x,x$adj_region,sd)->sd
table(x$adj_region)->n
t=m/sd
p=pt(-abs(t),df=n-1)*2
(!is.na(p) & m>m.cut[1] & p<alpha) -> sig; sig=as.numeric(sig) 
p=cbind(m,p,sig)
x=merge(x,p,by.x="adj_region",by.y="row.names")
x
}

##################

window.test=function(x=smo,alpha=0.05,m.cut=0)
{
p=pnorm(-abs(x$x)*2)
(!is.na(p) & x$x>m.cut & p<=alpha) -> sig; sig=as.numeric(sig) 
x=cbind(x,p,sig)
x
}

##############################################################

Adjust_Sig_Rigions=function(x,head=5,tail=5,ignore=5)
{
r=x$adj_region;o=x$order;s=x$sig
rs=table(r[s==1]);rs=rs[rs>ignore]
for (rsi in names(rs))
{
id=o[r==rsi];firstid=id[1];lastid=id[length(id)]
newlastid=lastid+tail; if (newlastid>length(o)){newlastid=length(o)}
s[lastid:newlastid]=1
}
x$adj_sig=s
x
}

##############################################################

Error_Of_Regions=function(sig,para)
{
#sig: significant regions, para: parameter matrix in simulation
sig0=numeric(length(sig))
sig0[para[3,"start"]:para[3,"end"]]=1
nh0=sum(sig0==0);nha=sum(sig0==1)
type1=sum(sig0==0 & sig==1)/nh0
power=sum(sig0==1 & sig==1)/nha
fdr=sum(sig0==0 & sig==1)/sum(sig==1)
sigid=c(1:length(sig))[sig==1]
t=para[3,"start"];t=sigid-t;b=t[order(abs(t))[1]]
t=para[3,"end"];t=sigid-t;a=t[order(abs(t))[1]]
err=cbind(type1,power,fdr,b,a)
err
}


##############################################################
# CN PLOTS ...
##############################################################

heatCN=function(x,co=gray(0:255/255),cut=NULL,xt=NULL,random=F)
{
ta=x
#raw data heatmap image
heatmap(ta,Colv=NA,scale="none")->ht
#if (random==FALSE) {ta=ta[ht$rowInd[1:357],]}
if (random==FALSE) {ta=ta[ht$rowInd,]}
if (random==TRUE) {ta=ta[sample(c(1:nrow(ta)),nrow(ta)),]}
ta[1,1]=0;ta[1,2]=6
if (length(cut)==2)
{
co=c("green","black","red")
bo=as.integer(min(ta)-5)
ta[ta<cut[1]]=bo-1;ta[ta>=cut[1] & ta<=cut[2]]=bo;ta[ta>cut[2]]=bo+1
ta=ta-bo
}
if (is.null(xt)){xt=paste(ncol(ta)," SNPs (Mbp)")}
xa=as.numeric(colnames(ta))/1000000
xa=seq(min(xa),max(xa),by=(max(xa)-min(xa))/length(xa))
image(x=xa,
y=c(1:nrow(ta)),z=t(ta),col=co,
main="Tumor/Normal Intensity Ratios",font=2,font.lab=2,
ylab=paste(nrow(ta)," Samples",sep=""),xlab=xt)
}
############################################################################
cn.plot=function(tt,type="raw",co=c(rgb(0,20:0/20,0),rgb(0:40/40,0,0)),
wsize=30,wn=800,grp=NULL,random=F)
{
temp=as.integer(ncol(tt)/wn+0.5)
if (temp<=0) {temp=1}
c=seq(1,ncol(tt),by=temp)
y=as.numeric(colnames(tt))/1000000
ta=tt[,c]
ta=2^ta*2
ta[ta>6]=6
xtext=paste(ncol(tt)," SNPs (Mbp)")

#-----CN heat map-----------------------
#png(paste(fi,"_raw.png",sep=""),height=600,width=800)
if (type=="raw")
{
if (is.null(grp)) heatCN(ta,co,xt=xtext,random=random)

if (!is.null(grp))
{
tb=ta[grp$id,]
#heatmap(tb,Colv=NA,scale="none")->ht
#b=tb[ht$rowInd,]
labr=c(1:nrow(tb))*NA;labr[1]=1;labr[length(labr)]=length(labr)
temp=seq(1,length(labr),by=as.integer(length(labr)/5));labr[temp]=temp
labc=c(1:ncol(tb))*NA;labc[1]=min(y[c]);labc[length(labc)]=max(y[c])
temp=seq(1,length(labc),by=as.integer(length(labc)/5))
labc[temp]=min(y[c])+(temp-1)*(max(y[c])-min(y[c]))/length(labc)
labc=substring(as.character(labc),1,5)
temp=unique(grp[,c("group","color")])
mtxt=character(0)
for (i in c(1:nrow(temp)))
{
mtxt=paste(mtxt,temp[i,"color"],":",temp[i,"group"],"   ",sep="")
}
heatmap(tb,Rowv=NA,Colv=NA,scale="none",col=co,
RowSideColors=as.character(grp$color),
labRow=labr,labCol=labc,cexRow=1,cexCol=1,
ylab="samples",xlab=xtext,main=mtxt)
}
}

#dev.off()
#png(paste(fi,"_2fd.png",sep=""),height=600,width=800)
#heatCN(ta,co,cut=c(0.5,2),xt=xtext)
#dev.off()

if (type=="CMDS")
{
if (is.null(grp))
{
#-----CMDS analysis -------------------- 

cmds(tt,wsize,10)->smo

#-----CMDS drawing
par(mfrow=c(2,2))
#---correlation matrix
image(x=y[c],y=y[c],z=cor(ta),
main="Correlation Matrix",font=2,font.lab=2,
ylab=xtext, xlab=xtext)
#---transformed r
plot(smo$pos,smo$x,type="l",font=2,font.lab=2,
main="Diagonal Tranformation",ylab="Transformed  r-value", xlab=xtext) #,ylim=c(-1,5))
#---segmentation
plot(smo$pos,smo$smth,type="l",font=2,font.lab=2,
main="Segmentation",ylab="Segmented  r-value", xlab=xtext) #,ylim=c(-1,5))
#---significant region(s)
plot(smo$pos,smo$sig,type="l",font=2,font.lab=2,yaxt="n",
main="Significant Region(s)",ylab="", xlab=xtext)
}


if (!is.null(grp))
{
par(mfrow=c(2,2))
for (gi in as.character(unique(grp[,"group"])))
{
idi=grp[,"id"][as.character(grp[,"group"])==gi]

Diagonalize_Regions(ta[idi,],w=wsize,ztrans=T)->temp
temp=cbind(c(1:length(temp)),temp)
Smooth_Regions(temp,nom=1,draw=F)-> smo
smo$adj_region=smo$region
Test_Regions(smo)->smo

#---correlation matrix
#image(x=y[c],y=y[c],z=cor(ta),
#main="Correlation Matrix",font=2,font.lab=2,
#ylab=xtext, xlab=xtext)
#---transformed r
#plot(y[c][1:nrow(smo)],smo$x,type="l",font=2,font.lab=2,
#main="Diagonal Tranformation",ylab="Transformed  r-value", xlab=xtext) #,ylim=c(-1,5))
#---segmentation
plot(y[c][1:nrow(smo)],smo$smth,type="l",font=2,font.lab=2,
main=gi,ylab="Segmented  r-value", xlab=xtext) #,ylim=c(-1,5))
#---significant region(s)
#plot(y[c][1:nrow(smo)],smo$sig,type="l",font=2,font.lab=2,yaxt="n",
#main="Significant Region(s)",ylab="", xlab=xtext)
}
}

}
}

# ----------------------- STAC simulation

to.stac=function(x,cut,DIR,ID,call.method="glad")
{
if (call.method=="glad")
{
tt=numeric(0)
for (i in c(1:nrow(x)))
{
temp=x[i,]
temp=cbind(c(1:length(temp)),temp)
Smooth_Regions(temp,norm=F)[,"smth"]->temp
tt=rbind(tt,temp)
}
}
if (call.method=="thresh") tt=x
if (call.method=="msa") tt=x
tt>=cut->tt
tt=tt*1
colnames(tt)=paste(c(1:ncol(tt)),"Mb",sep="")
rownames(tt)=paste("Exp",c(1:nrow(tt)),sep="")
invisible(tt)
if (!file.exists(DIR)) dir.create(DIR)
outf=paste(DIR,"/smiu",ID,".txt",sep="")
write.table(tt,quote=F,sep="\t",file=outf)
readLines(outf)->tt
tt[1]=paste("\t",tt[1],sep="")
write.table(tt,quote=F,file=outf,row.names=F,col.names=F)
}

#----------------------- STAC results

from.stac=function(DIR,PAR)
{
#DIR=outf;#PAR=param
dir(DIR)->fs
fs=paste(DIR,fs[grep("report",fs)],sep="/")
tt=numeric(0)
for (fi in fs)
{
read.table(fi,skip=1)$V4->sig0
for (alpha in c(10^-(1:10),5*10^-(1:10),(1:10)/10))
{
sig= (sig0<=alpha)*1
Error_Of_Regions(sig,PAR)->err; B=abs(err[1,"b"]);A=abs(err[1,"a"])
tt=rbind(tt,cbind(err,B,A,alpha))
}
}
by(tt,tt[,"alpha"],colMeans,na.rm=TRUE)->ttm
tt=numeric(0); for (ai in names(ttm)) tt=rbind(tt,ttm[[ai]]) 
tt
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

gvalues.test=function(x,sd.range=5,tail="two")
{
x=as.numeric(as.character(x))
sdx=abs(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
sd=(x-mean(x[sdx<sd.range],na.rm=TRUE))/sd(x[sdx<sd.range],na.rm=TRUE)
if (tail=="two") p=2*pnorm(-abs(sd))
if (tail=="right") p=pnorm(-sd)
if (tail=="left") p=pnorm(sd)
tt=cbind(sd,p)
invisible(tt)
}
#######################################

cmds.v1=function(x,wsize=10,wstep=1,ztrans=T,sd.range=5)
{
Diagonalize_Regions(x,wsize,wstep,ztrans)->z
gvalues.test(x=z[,"z"],sd.range=sd.range,tail="right")->zt; colnames(zt)=paste("z",colnames(zt),sep=".")
gvalues.test(x=z[,"m"],sd.range=sd.range,tail="two")->mt; colnames(mt)=paste("m",colnames(mt),sep=".")
z=cbind(z,mt,zt)
invisible(z)
}

####################################################### Diagonalize_Regions()

Diagonalize_Regions=function(x,w=30,wby=1,ztrans=TRUE)
{
n=ncol(x)
z=seq(1,(n-w),by=wby); m=z
window=c(1:length(z));start=z;mid=z+as.integer(w/2); end=z+w
for (i in window)
{
xi=x[,z[i]:(z[i]+w)]
m[i]= mean(xi,na.rm=TRUE)
ci=cor(xi,use="pairwise.complete.obs")
for (ii in c(1:dim(ci)[1])){ci[ii,ii]=NA}
if (ztrans==TRUE)
{
ci=0.5*log((1+ci)/(1-ci))
ci[abs(ci)==Inf]=NA
ci=ci*sqrt(nrow(x)-3)
}
z[i]= mean(ci,na.rm=TRUE)
}
z=cbind(window,start,mid,end,m,z)
}

################################################ cmds.focal.test()

cmds.focal.test=function(
data.dir="your_data_dir",
wsize=30,
wstep=1,
analysis.ID=NA,
chr.colname="chromosome",
pos.colname="position", 
plot.dir="focal_plot",
result.dir="focal") 

{
dir.create(plot.dir,recursive=T)
dir.create(result.dir,recursive=T)
fs=dir(data.dir); if (!is.na(analysis.ID)) fs=fs[as.numeric(analysis.ID)]

for (fi in fs)
{
in.file=paste(data.dir,fi,sep="/")
out.file=paste(result.dir,"/",fi,".test",sep="")
out.image=paste(plot.dir,"/",fi,".png",sep="")
cat("file is ");
cat(in.file)
read.table(in.file,header=T)->x
chromosome=unique(as.character(x[,chr.colname]))
pos=x[,pos.colname]
x=t(x[,-(1:2)])
cmds.v1(x,wsize=wsize,wstep=wstep)->smo
smo=as.numeric.frame(smo,5:10,5)
smo[,"start"]=pos[smo[,"start"]]
smo[,"mid"]=pos[smo[,"mid"]]
smo[,"end"]=pos[smo[,"end"]]
smo=cbind(chromosome,smo)
write.table(smo,file=out.file,row.names=F,quote=F,sep="\t")
bitmap(out.image,width=1024,height=768,units="px");cmds.test.plot(smo);dev.off()
}

}

################################## cmds.test.plot()

cmds.test.plot=function(smo,yscale=NULL)
{
par(mfcol=c(2,2))
plot(smo[,"mid"],smo[,"m.sd"],type="l",lwd=2,main="Mean CN of All Samples",
xlab="Position",ylab="m.sd")
abline(0,0,col="red")

plot(smo[,"mid"],-log10(smo[,"m.p"]),type="l",lwd=2,main="Test of Mean (H0:m.sd=0)",
xlab="Position",ylab="-log(P)")

plot(smo[,"mid"],smo[,"z.sd"],type="l",lwd=2,main="RCNA Score",
xlab="Position", ylab="z.sd")
abline(0,0,col="red")

plot(smo[,"mid"],-log10(smo[,"z.p"]),type="l",lwd=2,main="CMDS Test (H0:z.sd=0)",
xlab="Position",ylab="-log(P)",ylim=yscale)
#abline(2,0,lty=2,col="red",lwd=2)
}
#######################################################################
whole_genome_test_plots=function(all_chr_data_file)
{
    x=read.table(all_chr_data_file,header=T);
    x$z.p = -log10(x$z.p);
    x$m.p = -log10(x$m.p);

    #plot mean test
    img_filename=paste(all_chr_data_file,".meantest",sep="");
    plotgenome(x,y="m.p",yscale=c(0,6),pos.col="mid",img=img_filename,chr.col="chromosome",suffix="m.p",cutline=NULL);


    #plot cmds test
    img_filename=paste(all_chr_data_file,".cmdstest",sep="");
    plotgenome(x,y="z.p",yscale=c(0,6),pos.col="mid",img=img_filename,chr.col="chromosome",suffix="cmds.p",cutline=NULL);
}
#######################################################################

whole_multi_genome_test_plots=function(all_chr_data_files,outfile,chrom=NULL)
{
    data=NULL;
    #all_chr_data_files is a comma-delimited list of files to plot on the same graphs
    chars=strsplit(as.character(all_chr_data_files),"");
    files=strsplit(all_chr_data_files,",");

    for(file in c(1:length(files[[1]])))
    {
        x=read.table(files[[1]][file],header=T);
        if (is.null(data)) {
            data=cbind(as.character(x$chromosome),as.numeric(x$mid));
            colnames(data)=c("CHR","POS");
        }
        z = -log10(x$z.p);
        data=cbind(data,z);
        colname=paste("file",file,sep="");
        colnames(data)[file+2]=colname;
    }
    print(tail(data))
    #stop();
    save(data,file="testdata_cmds_plotmulti");

    #plot cmds test
    #img_filename=paste(all_chr_data_file,".cmdstest",sep="");
    plot_multi_genome(as.matrix(data),y=colnames(data)[3:length(colnames(data))],num_data_cols=length(files[[1]]),yscale=c(0,6),pos.col="POS",img=outfile,chr.col="CHR",suffix="",cutline=NULL,chrom=chrom);
}



#####################################################################################
plot_multi_genome = function (tt,y="p",num_data_cols=1,cutoff=NULL,cutline=2,img=NULL,yscale=NULL,draw=TRUE,ltype="p",
chrom=NULL,mbp=NULL,chr.col="chromosome",pos.col="position",tombp=T,suffix="")
{
#colors: blue, red, green,orange,pink,purple,black,yellow,cyan,brown
colors_used=1:num_data_cols;
colnames(tt)[grep(chr.col,colnames(tt))]="chr"
colnames(tt)[grep(pos.col,colnames(tt))]="mbp"
tt=as.data.frame(tt);
tt$mbp=as.numeric(as.character(tt$mbp))
for (coli in y)
{
    tt[,coli] = as.numeric(as.character(tt[,coli]));
}

if (tombp) tt$mbp=tt$mbp/1000000
if (is.null(chrom)) {chrom=levels(as.factor(tt[,"chr"]));chrom=intersect(c(1:24,"X","Y","x","y"),chrom)}
print(chrom)
if (length(mbp)==2) {tt=tt[(tt[,"mbp"]>=mbp[1] & tt[,"mbp"]<=mbp[2]),]}
if (!is.null(cutoff)) {tt[,y][tt[,y]<cutoff]=NA}
if (draw==TRUE)
{
#if (length(img)>0) {png(paste(img,"png",sep="."),1200,800)}
if (length(img)>0) {bitmap(paste(img,"png",sep="."),width=1200,height=800,units="px")}

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
cat(paste("in for chr in chrom",chr))
tt[,"chr"]==chr->ch
tl=paste("Chrom.",chr,suffix,sep=" ")
matplot(tt[ch,"mbp"],tt[ch,y],pch=".",main=tl,xlab="Position (Mbp)",ylab="-log(P)",ylim=yscale,type=ltype,col=colors_used)
for (ct in cutline)
{abline(ct,0,col="red")}
}#chr
if (length(img)>0) {dev.off()}
}#draw
tt=tt[!is.na(tt[,y]),]
invisible(tt)
}

#####################################################################
plot_chr_vertical=function(data_file)
{
#before plotting, you need to save the whole genome merged data into one file
#e.g. cmdsAll.rdata,  using data object name "cmds"

output_filename = paste(data_file,".vertchr.png",sep="");
read.table(data_file,header=TRUE)->cmds;
#(load("/gscmnt/sata423/info/medseq/llin/aml/copy_number_analysis/all_Jackie/CMDS/cmdsFocalAll.rdata"));
#cmds=tt;
#head(cmds)


#png("CMDS150_AllinOne.png");
png(output_filename);
#y=cmds[,c("primary","met")];y.lab="Log2Ratio";tl="Primary(blue) vs Metastasis(red)";lim=c(-4,4);

#co=c("blue","red");
y=-log10(cmds[,"z.p"]);y.lab="-log(P)";tl="CMDS test";lim=c(-20,20);coamp="red";codel="blue";
yamp=y;

ydel=log10(cmds[,"z.p"]);

yamp[yamp>lim[2]]=lim[2];
ydel[ydel<lim[1]]=lim[1];
yamp=yamp-3;
ydel=ydel+3;

yamp[cmds$m<0]="NA";
ydel[cmds$m>0]="NA";

yamp[yamp<0]="NA";
ydel[ydel>0]="NA";


#y=log2(cmds[,"m"]/2);y.lab="Log2Ratio";tl="Mean";lim=c(-0.5,0.5);co="black";
##### start .....

chr=as.character(cmds[,"chromosome"])
mbp=cmds[,"mid"]/1000000
len=0
chr.names=intersect(c(1:24),unique(chr))
aclen=c(1:(length(chr.names)+1))
names(aclen)=c(0,chr.names)
for (ci in c(1:length(chr.names)))
{
len=c(len,max(mbp[chr==chr.names[ci]]))
aclen[ci]=sum(len[1:ci])
mbp[chr==chr.names[ci]]= mbp[chr==chr.names[ci]]+ aclen[ci]
}
aclen[length(aclen)]=sum(len[1:length(aclen)])
#y=y[order(mbp)];
#y[y>lim[2]]=lim[2];y[y<lim[1]]=lim[1]
mbp=mbp[order(mbp)]

#### draw ......


#matplot(mbp,y,type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=co,axes=FALSE,add=0)
#matplot(mbp,yamp,type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=coamp,axes=FALSE,add=0)
#matplot(mbp,ydel,type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=codel,axes=FALSE,add=0)
##matplot(mbp,cbind(yamp,ydel),type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=c("red","blue"),axes=FALSE,add=0)
matplot(cbind(yamp,ydel),mbp,type="l",main=tl,xlab=y.lab,ylab="",xlim=lim,col=c("red","blue"),axes=FALSE,add=0)
#lines(c(0,aclen[length(aclen)]),c(0,0))
#for (i in seq(1,length(chr.names),by=2))
#{
#xi=c(aclen[i],aclen[i+1],aclen[i+1],aclen[i])
#yi=c(lim[1],lim[1],lim[2],lim[2])
#polygon(xi,yi,border=NA,col=grey(0.9))
#}
lim=c(-17,17);
for (i in seq(1,length(chr.names),by=2))
{
yi=c(aclen[i],aclen[i+1],aclen[i+1],aclen[i])
xi=c(lim[1],lim[1],lim[2],lim[2])
polygon(xi,yi,border=NA,col=grey(0.9))
}


#matplot(mbp,yamp,type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=coamp,axes=FALSE,add=1)
#matplot(mbp,ydel,type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=codel,axes=FALSE,add=1)
##matplot(mbp,cbind(yamp,ydel),type="l",main=tl,xlab="Chromosomes",ylab=y.lab,ylim=lim,col=c("red","blue"),axes=FALSE,add=1)
matplot(cbind(yamp,ydel),mbp,type="l",main=tl,ylab="",xlab=y.lab,xlim=lim,col=c("red","blue"),axes=FALSE,add=1)
#ylbs=seq(lim[1],lim[2],by=(lim[2]-lim[1])/4)
#ylbs=c(20,10,0,10,20)
#xlbs=c(20,10,0,10,20)
xlbs=c(20,10,5,3,5,10,20)
axis(side=1,at=c(17,7,2,0,-2,-7,-17),labels=xlbs)
#axis(side=2,at=ylbs,labels=ylbs)
#axis(side=2,at=c(20,10,0,-10,-20),labels=ylbs)
#axis(side=1,at=c(20,10,0,-10,-20),labels=xlbs)
#text(x=((aclen[2:length(aclen)]-aclen[1:length(chr.names)])/2+aclen[1:length(chr.names)]),y=rep((lim[1]-1),times=length(chr.names)),labels=chr.names,cex=0.5)
text(y=((aclen[2:length(aclen)]-aclen[1:length(chr.names)])/2+aclen[1:length(chr.names)]),x=rep((lim[1]-1),times=length(chr.names)),labels=chr.names,cex=0.5)

dev.off()
}
################################### plotgenome()

plotgenome = function (tt, y="p",cutoff=NULL,cutline=2,img=NULL,yscale=NULL,draw=TRUE,ltype="p",
chrom=NULL,mbp=NULL,chr.col="chromosome",pos.col="position",tombp=T,suffix="")
{
colnames(tt)[grep(chr.col,colnames(tt))]="chr"
colnames(tt)[grep(pos.col,colnames(tt))]="mbp"
if (tombp) tt$mbp=tt$mbp/1000000
if (is.null(chrom)) {chrom=levels(as.factor(tt[,"chr"]));chrom=intersect(c(1:24,"X","Y","x","y"),chrom)}
if (length(mbp)==2) {tt=tt[(tt[,"mbp"]>=mbp[1] & tt[,"mbp"]<=mbp[2]),]}
if (!is.null(cutoff)) {tt[,y][tt[,y]<cutoff]=NA}
if (draw==TRUE)
{
if (length(img)>0) {bitmap(paste(img,"png",sep="."),width=1200,height=800,units="px")}
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
tl=paste("Chrom.",chr,suffix,sep=" ")
matplot(tt[ch,"mbp"],tt[ch,y],pch=".",main=tl,xlab="Position (Mbp)",ylab=y,ylim=yscale,type=ltype)
for (ct in cutline) {abline(ct,0,col="red")}
}#chr
if (length(img)>0) {dev.off()}
}#draw

tt=tt[!is.na(tt[,y]),]
invisible(tt)
}

##################################################### Region_calls()

Region_calls = function(datafile,chr,start,end,permun=NA,output_dir=".")

{

start=as.integer(start);
end=as.integer(end);
#end=as.numeric(end);
permun=as.numeric(permun);

read.table(datafile,header=T,sep="\t",comment.char="")->x
head(x)
dim(x)

z=x[as.character(x[,1])==chr & x[,2]>=start & x[,2]<=end,]
head(x)
dim(x)

samples=colnames(z)[-(1:2)]
cn=colMeans(z[,-(1:2)],na.rm=T)
mean=mean(cn,na.rm=T)

#sd=cn-mean/sd(cn,na.rm=T)
sd=(cn-mean)/sd(cn,na.rm=T)

p=pt(-abs(sd),ncol(z)-3)
fdr=p.adjust(p,method="BH")
pp=cbind(samples,mean,cn,sd,p,fdr)

if (!is.na(permun))
{
samplen=nrow(z)
totaln=c(1:nrow(x))
x=x[,-(1:2)]
pa=cn*0; pd=pa
for (i in 1:permun)
{
z=x[sample(totaln,samplen),]
cni=colMeans(z)
pa=pa + as.numeric(cni>=cn)
pd=pd + as.numeric(cni<=cn)
}

#p.permu=pa
#p.permu[pd<pa]=pd[pd<pa]
#p.permu=p.permu/permun

p.permu=pa*NA
p.permu[sd<=0]=pd[sd<=0]
p.permu[sd>0]=pa[sd>0]

fdr.permu=p.adjust(p.permu,method="BH")
pp=cbind(pp,p.permu,fdr.permu)
}

pp=pp[order(cn),]
outfile_chr = paste(output_dir,chr,sep="/");
write.csv(pp,file=paste(outfile_chr,start,end,"call.csv",sep="_"),quote=F,row.names=F)

}

####################################################################################

#Copy number data segmentation using CBS algorithm in R BioConductor package DNAcopy

#example data
#cn.file="/gscmnt/sata181/info/medseq/llin/BC/SNP/Infinium1MOmni_EllisBreastCancer_20100125_analysis/CN/BRC9T.log2"
#map.file="/gscmnt/sata181/info/medseq/llin/BC/SNP/Infinium1MOmni_EllisBreastCancer_20100125_analysis/map.csv"

cbs=function(cn.file,map.file,output.dir=".")
{
temp=strsplit(cn.file,split="/")[[1]]
temp=temp[length(temp)]
output.file=paste(output.dir,"/",temp,"_cbs",sep="")
library(DNAcopy)
read.table(cn.file,header=T,nrows=-1)->cn
read.table(map.file,header=T,nrow=-1)->map
cno <- CNA(as.numeric(cn[,1]),map[,1],map[,2],data.type="logratio",sampleid=names(cn))
seg <- segment(smooth.CNA(cno), undo.splits = "sdundo", undo.SD=2)
seg = seg$output[,-1]
#plot(seg,plot.type="s")
write.table(seg, file=output.file, row.names=F, quote=F)
invisible(seg)
}

