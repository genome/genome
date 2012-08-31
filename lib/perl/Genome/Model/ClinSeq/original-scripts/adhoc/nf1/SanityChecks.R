#Targeted hypothesis testing (sanity checks)

#Load libraries
library(RColorBrewer)
library(genefilter)
library(scatterplot3d)
library(gcrma)
library(preprocessCore)

#Simple examination of variability between replicates vs across conditions
setwd("/Users/mgriffit/Documents/WASHU/Projects/Josh_Rubin")
dir()

#Identify the classes of each data column
male_female=c("M","M","M","F","F","F","M","M","M","F","F","F")
table(male_female)
wt_null=c("Wt","Wt","Wt","Wt","Wt","Wt","Null","Null","Null","Null","Null","Null")
table(wt_null)
full_class=c("Wt_M","Wt_M","Wt_M","Wt_F","Wt_F","Wt_F","Null_M","Null_M","Null_M","Null_F","Null_F","Null_F")
table(full_class)
new_order=c("Wt_M_1", "Wt_M_2", "Wt_M_3", "Wt_F_1", "Wt_F_2", "Wt_F_3", "Null_M_1", "Null_M_2", "Null_M_3", "Null_F_1", "Null_F_2", "Null_F_3")

cols=unclass(factor(male_female))[1:length(male_female)]
colors=brewer.pal(4,"Set1")
male_female_cols=colors[cols]

cols=unclass(factor(wt_null))[1:length(wt_null)]
colors=brewer.pal(4,"Set1")
wt_null_cols=colors[cols]

cols=unclass(factor(full_class))[1:length(full_class)]
colors=brewer.pal(4,"Set1")
class_cols=colors[cols]


#Load the data Y chromosome data
y_data=read.table("Y_genes.txt", sep="\t", header=TRUE, as.is=1:2)
names(y_data)=c("tracking_id", "gene_names", "Wt_M_1", "Wt_F_1", "Null_M_1", "Null_F_1", "Null_F_2", "Null_F_3", "Wt_M_3", "Wt_F_2", "Wt_F_3", "Null_M_3", "Null_M_2", "Wt_M_2")

y_fpkm=y_data[,new_order]

#Set 0's to NA to avoid log2 of 0 instead
y_fpkm2=y_fpkm
y_fpkm2[y_fpkm==0]=NA
boxplot(log2(y_fpkm2), col=class_cols, ylab="Log2 FPKM, zeros removed", main="Distribution of *ChrY* FPKMs for all samples", las=2)

#Calculate a correlation matrix for all genes after removing 0's
r=cor(y_fpkm, use="pairwise.complete.obs", method="spearman")

#set na's to 0
r[is.na(r)]=0
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all *ChrY* genes)")
text(mds$points[,1], mds$points[,2], male_female, col=male_female_cols)




#Load the data X chromosome data
x_data=read.table("X_genes.txt", sep="\t", header=TRUE, as.is=1:2)
names(x_data)=c("tracking_id", "gene_names", "Wt_M_1", "Wt_F_1", "Null_M_1", "Null_F_1", "Null_F_2", "Null_F_3", "Wt_M_3", "Wt_F_2", "Wt_F_3", "Null_M_3", "Null_M_2", "Wt_M_2")
x_fpkm=x_data[,new_order]

#Set 0's to NA to avoid log2 of 0 instead
x_fpkm2=x_fpkm
x_fpkm2[x_fpkm==0]=NA
boxplot(log2(x_fpkm2), col=class_cols, ylab="Log2 FPKM, zeros removed", main="Distribution of *ChrX* FPKMs for all samples", las=2)

#Calculate a correlation matrix for all genes after removing 0's
r=cor(x_fpkm, use="pairwise.complete.obs", method="spearman")

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all *ChrX* genes)")
text(mds$points[,1], mds$points[,2], male_female, col=male_female_cols)

#Redo with a coefficient of variation filter instead of straight variance
min_cov=1
mmfilt=function(minv=min_variance){
	function(x){
		cov=sd(x, na.rm=TRUE)/mean(x,na.rm=TRUE)
		(cov > min_cov)
	}	
}
mmfun = mmfilt()
ffun = filterfun(mmfun)
i = genefilter(x_fpkm, ffun)
sum(i)
r=cor(x_fpkm[i,], use="pairwise.complete.obs", method="spearman")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
min_x=min(mds$points[,1])-abs((min(mds$points[,1])*zz))
max_x=max(mds$points[,1])+abs((max(mds$points[,1])*zz))
min_y=min(mds$points[,2])-abs((min(mds$points[,2])*zz))
max_y=max(mds$points[,2])+abs((max(mds$points[,2])*zz))
par(mfrow=c(2,2))
main=paste("(cov > ", min_cov, "; n genes = ", sum(i), ")", sep="")
plot(mds$points, type="n", xlab="", ylab="", main=paste("Samples", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], new_order, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Classes", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], full_class, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("WT/Null", main)); text(mds$points[,1], mds$points[,2], wt_null, col=wt_null_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Male/Female", main)); text(mds$points[,1], mds$points[,2], male_female, col=male_female_cols)











