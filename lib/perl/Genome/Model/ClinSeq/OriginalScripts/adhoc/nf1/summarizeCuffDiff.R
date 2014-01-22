library(RColorBrewer)
library(genefilter)
library(scatterplot3d)
library(gcrma)
library(preprocessCore)

#Initially used the following data as input:
#/gscmnt/sata132/techd/solexa/jwalker/RNAseq/Mouse/Josh_Rubin/Josh_Rubin_ab_initio/cuffdiff/genes.fpkm_tracking.simple

#Try repeating with:
#/gscmnt/sata132/techd/solexa/jwalker/RNAseq/Mouse/Josh_Rubin/Josh_Rubin_ab_initio/cuffdiff_wGene/genes.fpkm_tracking.simple


#Simple examination of variability between replicates vs across conditions
setwd("/Users/mgriffit/Documents/WASHU/Projects/Josh_Rubin")
dir()

#Load the data
data=read.table("genes.fpkm_tracking.simple", sep="\t", header=TRUE, as.is=1:2)

#Assign more readable names to the columns of the data frame
names(data)=c("tracking_id", "gene_names", "Wt_M_1", "Wt_F_1", "Null_M_1", "Null_F_1", "Null_F_2", "Null_F_3", "Wt_M_3", "Wt_F_2", "Wt_F_3", "Null_M_3", "Null_M_2", "Wt_M_2")


#Set a more human readable order for the columns and grab a matrix with just the data values
new_order=c("Wt_M_1", "Wt_M_2", "Wt_M_3", "Wt_F_1", "Wt_F_2", "Wt_F_3", "Null_M_1", "Null_M_2", "Null_M_3", "Null_F_1", "Null_F_2", "Null_F_3")
fpkm=data[,new_order]

#Identify the classes of each data column
male_female=c("M","M","M","F","F","F","M","M","M","F","F","F")
table(male_female)
wt_null=c("Wt","Wt","Wt","Wt","Wt","Wt","Null","Null","Null","Null","Null","Null")
table(wt_null)
full_class=c("Wt_M","Wt_M","Wt_M","Wt_F","Wt_F","Wt_F","Null_M","Null_M","Null_M","Null_F","Null_F","Null_F")
table(full_class)

#Assign colors to each column based on their classes
display.brewer.all()

cols=unclass(factor(male_female))[1:length(male_female)]
colors=brewer.pal(4,"Set1")
male_female_cols=colors[cols]

cols=unclass(factor(wt_null))[1:length(wt_null)]
colors=brewer.pal(4,"Set1")
wt_null_cols=colors[cols]

cols=unclass(factor(full_class))[1:length(full_class)]
colors=brewer.pal(4,"Set1")
class_cols=colors[cols]

#Look at the distribution of FPKM values before getting carried away
#Convert to log2 scale to deal with outliers and add 1 to prevent log2(0)
boxplot(log2(fpkm+1), col=class_cols, ylab="Log2 FPKM+1", main="Distribution of FPKMs for all samples", las=2)

#Set 0's to NA to avoid log2 of 0 instead
fpkm2=fpkm
fpkm2[fpkm==0]=NA
boxplot(log2(fpkm2), col=class_cols, ylab="Log2 FPKM, zeros removed", main="Distribution of FPKMs for all samples", las=2)

#Calculate a correlation matrix for all genes after removing 0's
r=cor(fpkm2, use="pairwise.complete.obs", method="spearman")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
plot(mds$points, type="n", xlab="", ylab="", xlim=c(-0.04, 0.05), ylim=c(-0.03,0.02), main="MDS distance plot (all non-zero genes)")
text(mds$points[,1], mds$points[,2], full_class, col=class_cols)

#plot space expansion factor
zz=0.2
#Now repeat the same analysis but identify those genes that are actually varying across the libraries
#Very simply, apply a variance filter
min_variance=10
mmfilt=function(minv=min_variance){
	function(x){
		v = var(x)
		(v > min_variance)
	}	
}
mmfun = mmfilt()
ffun = filterfun(mmfun)
i = genefilter(fpkm, ffun)
sum(i)
r=cor(fpkm[i,], use="pairwise.complete.obs", method="spearman")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
min_x=min(mds$points[,1])-abs((min(mds$points[,1])*zz))
max_x=max(mds$points[,1])+abs((max(mds$points[,1])*zz))
min_y=min(mds$points[,2])-abs((min(mds$points[,2])*zz))
max_y=max(mds$points[,2])+abs((max(mds$points[,2])*zz))
par(mfrow=c(2,2))
main=paste("(var > ", min_variance, "; n genes = ", sum(i), ")", sep="")
plot(mds$points, type="n", xlab="", ylab="", main=paste("Samples", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], new_order, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Classes", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], full_class, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("WT/Null", main)); text(mds$points[,1], mds$points[,2], wt_null, col=wt_null_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Male/Female", main)); text(mds$points[,1], mds$points[,2], male_female, col=male_female_cols)


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
i = genefilter(fpkm, ffun)
sum(i)
r=cor(fpkm[i,], use="pairwise.complete.obs", method="spearman")
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

#Now try it in 3d space
mds=cmdscale(d, k=3, eig=TRUE)
scatterplot3d(x=mds$points[,1], y=mds$points[,2], z=mds$points[,3], xlab="", ylab="", zlab="", color=class_cols, pch=16)
scatterplot3d(x=mds$points[,1], y=mds$points[,2], z=mds$points[,3], xlab="", ylab="", zlab="", color=wt_null_cols, pch=16)
scatterplot3d(x=mds$points[,1], y=mds$points[,2], z=mds$points[,3], xlab="", ylab="", zlab="", color=male_female_cols, pch=16)

#Try a min expression level?
min_exp=1
mmfilt=function(minv=min_exp){
	function(x){
		s = sum(x)
		(s > min_exp)
	}	
}
mmfun = mmfilt()
ffun = filterfun(mmfun)
i = genefilter(fpkm, ffun)

sum(i)
r=cor(fpkm[i,], use="pairwise.complete.obs", method="spearman")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
min_x=min(mds$points[,1])-abs((min(mds$points[,1])*zz))
max_x=max(mds$points[,1])+abs((max(mds$points[,1])*zz))
min_y=min(mds$points[,2])-abs((min(mds$points[,2])*zz))
max_y=max(mds$points[,2])+abs((max(mds$points[,2])*zz))
par(mfrow=c(2,2))
main=paste("(exp > ", min_exp, "; n genes = ", sum(i), ")", sep="")
plot(mds$points, type="n", xlab="", ylab="", main=paste("Samples", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], new_order, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Classes", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], full_class, col=class_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("WT/Null", main)); text(mds$points[,1], mds$points[,2], wt_null, col=wt_null_cols)
plot(mds$points, type="n", xlab="", ylab="", main=paste("Male/Female", main)); text(mds$points[,1], mds$points[,2], male_female, col=male_female_cols)

#normalize the data to allow head-to-head correlations
names=names(fpkm)
x = as.matrix(fpkm)
fpkm_norm = normalize.quantiles(x)
fpkm_norm = as.data.frame(fpkm_norm)
dimnames(fpkm_norm)[[2]] = names

#Look at the distribution of FPKM values before getting carried away
#Convert to log2 scale to deal with outliers and add 1 to prevent log2(0)
par(mfrow=c(2,2))
boxplot(log2(fpkm+1), col=class_cols, ylab="Log2 FPKM+1", xlab="", main="FPKMs (Transformed only)", las=2)
boxplot(log2(fpkm_norm+1), col=class_cols, ylab="Log2 FPKM_Norm+1", xlab="", main="FPKMs (Transformed & Normalized)", las=2)

#Set 0's to NA to avoid log2 of 0 instead
fpkm2=fpkm
fpkm2[fpkm==0]=NA
fpkm_norm2=fpkm_norm
fpkm_norm2[fpkm_norm==0]=NA
boxplot(log2(fpkm2), col=class_cols, ylab="Log2 FPKM, zeros removed", xlab="", main="FPKMs (Transformed only)", las=2)
boxplot(log2(fpkm_norm2), col=class_cols, ylab="Log2 FPKM_Norm, zeros removed", xlab="", main="FPKMs (Transformed & Normalized)", las=2)

#Pick a few pairwise examples to plot against each other
legend_text = c("No change","Fold change > 2")
legend_lwd=2
legend_cols=c("red", "orange")
legend_lty=2

plot(x=log2(fpkm_norm2[,"Wt_M_2"]), y=log2(fpkm_norm2[,"Wt_M_3"]), xlim=c(-6,15), ylim=c(-6,15), main="Correlation between FPKMs for a pair of libraries", xlab="WT, Male, Lot2", ylab="WT, Male, Lot3")
abline(0,1, col="red")
abline(1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
abline(-1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
legend("topleft", legend_text, lwd=legend_lwd, lty=legend_lty, col=legend_cols)

plot(x=log2(fpkm_norm2[,"Wt_M_3"]), y=log2(fpkm_norm2[,"Null_M_3"]), xlim=c(-6,15), ylim=c(-6,15), main="Correlation between FPKMs for a pair of libraries", xlab="WT, Male, Lot3", ylab="Null, Male, Lot3")
abline(0,1, col="red")
abline(1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
abline(-1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
legend("topleft", legend_text, lwd=legend_lwd, lty=legend_lty, col=legend_cols)

#Add an arbitrary number to avoid log2 of 0
min_nonzero=0.0001

####################################################
#Plot all within lot comparisons.  Include 0s.
plotComp=function(indata, s1, s1_name, s2, s2_name, main_name){
	#Get min non-zero value from the two vectors
	z=indata
	z[indata==0]=NA
	min1=min(z[,s1], na.rm=TRUE)	
	min2=min(z[,s2], na.rm=TRUE)	
	min_nonzero=min(c(min1,min2))
	legend_text = c("No change","Fold change > 2", "Fold change > 4")
	legend_lwd=2
	legend_cols=c("black", "orange", "red")
	legend_lty=2
	plot(x=log2(indata[,s1]+min_nonzero), y=log2(indata[,s2]+min_nonzero), main=main_name, 
		 xlab=s1_name, ylab=s2_name, col="grey20")
	abline(0,1, col="black", lwd=2)
	abline(1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
	abline(-1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
	abline(2,1, col="red", lty=legend_lty, lwd=legend_lwd)
	abline(-2,1, col="red", lty=legend_lty, lwd=legend_lwd)
	legend("topleft", legend_text, lwd=legend_lwd, lty=c(1,legend_lty,legend_lty), col=legend_cols)
}


ttest=function(indata){
    obj<-try(t.test(x=indata[set1], y=indata[set2]), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

plotComp(fpkm, "Wt_M_1", "FPKM (WT, Male, Lot1)", "Wt_F_1", "FPKM (WT, Female, Lot1)", "Correlation between FPKMs, WT Male vs. Female, Lot1")
plotComp(fpkm, "Wt_M_2", "FPKM (WT, Male, Lot2)", "Wt_F_2", "FPKM (WT, Female, Lot2)", "Correlation between FPKMs, WT Male vs. Female, Lot2")
plotComp(fpkm, "Wt_M_3", "FPKM (WT, Male, Lot3)", "Wt_F_3", "FPKM (WT, Female, Lot3)", "Correlation between FPKMs, WT Male vs. Female, Lot3")
plotComp(fpkm, "Null_M_1", "FPKM (Null, Male, Lot1)", "Null_F_1", "FPKM (Null, Female, Lot1)", "Correlation between FPKMs, NULL Male vs. Female, Lot1")
plotComp(fpkm, "Null_M_2", "FPKM (Null, Male, Lot2)", "Null_F_2", "FPKM (Null, Female, Lot2)", "Correlation between FPKMs, NULL Male vs. Female, Lot2")
plotComp(fpkm, "Null_M_3", "FPKM (Null, Male, Lot3)", "Null_F_3", "FPKM (Null, Female, Lot3)", "Correlation between FPKMs, NULL Male vs. Female, Lot3")

plotComp(fpkm, "Wt_M_1", "FPKM (WT, Male, Lot1)", "Null_M_1", "FPKM (Null, Male, Lot1)", "Correlation between FPKMs, Male WT vs. NULL, Lot1")
plotComp(fpkm, "Wt_M_2", "FPKM (WT, Male, Lot2)", "Null_M_2", "FPKM (Null, Male, Lot2)", "Correlation between FPKMs, Male WT vs. NULL, Lot2")
plotComp(fpkm, "Wt_M_3", "FPKM (WT, Male, Lot3)", "Null_M_3", "FPKM (Null, Male, Lot3)", "Correlation between FPKMs, Male WT vs. NULL, Lot3")
plotComp(fpkm, "Wt_F_1", "FPKM (WT, Female, Lot1)", "Null_F_1", "FPKM (Null, Female, Lot1)", "Correlation between FPKMs, Female WT vs. NULL, Lot1")
plotComp(fpkm, "Wt_F_2", "FPKM (WT, Female, Lot2)", "Null_F_2", "FPKM (Null, Female, Lot2)", "Correlation between FPKMs, Female WT vs. NULL, Lot2")
plotComp(fpkm, "Wt_F_3", "FPKM (WT, Female, Lot3)", "Null_F_3", "FPKM (Null, Female, Lot3)", "Correlation between FPKMs, Female WT vs. NULL, Lot3")

#Calculate the mean FPKM across replicates so that these graphs can be summarized further
#Use the normalized data for these comparisons
fpkm_norm[,"Wt_M"] = apply(fpkm_norm[,c("Wt_M_1","Wt_M_2","Wt_M_3")], 1, mean)
fpkm_norm[,"Wt_F"] = apply(fpkm_norm[,c("Wt_F_1","Wt_F_2","Wt_F_3")], 1, mean)
fpkm_norm[,"Null_M"] = apply(fpkm_norm[,c("Null_M_1","Null_M_2","Null_M_3")], 1, mean)
fpkm_norm[,"Null_F"] = apply(fpkm_norm[,c("Null_F_1","Null_F_2","Null_F_3")], 1, mean)
fpkm_norm[,"M"] = apply(fpkm_norm[,c("Wt_M_1","Wt_M_2","Wt_M_3","Null_M_1","Null_M_2","Null_M_3")], 1, mean)
fpkm_norm[,"F"] = apply(fpkm_norm[,c("Wt_F_1","Wt_F_2","Wt_F_3","Null_F_1","Null_F_2","Null_F_3")], 1, mean)
fpkm_norm[,"Wt"] = apply(fpkm_norm[,c("Wt_M_1","Wt_M_2","Wt_M_3","Wt_F_1","Wt_F_2","Wt_F_3")], 1, mean)
fpkm_norm[,"Null"] = apply(fpkm_norm[,c("Null_M_1","Null_M_2","Null_M_3","Null_F_1","Null_F_2","Null_F_3")], 1, mean)

#Calculate p-values and log2 foldchanges for all comparisons.  Add an arbitrary number to stabilize variance, and reduce the effect of large differences in very lowly expressed genes
stab_var = 0.1
set1=c("Wt_M_1","Wt_M_2","Wt_M_3"); set2=c("Wt_F_1","Wt_F_2","Wt_F_3")
fpkm_norm[,"Wt_M_vs_F_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"Wt_M_vs_F_diff"] = log2(fpkm_norm[,"Wt_M"] + stab_var) - log2(fpkm_norm[,"Wt_F"] + stab_var)

set1=c("Null_M_1","Null_M_2","Null_M_3"); set2=c("Null_F_1","Null_F_2","Null_F_3")
fpkm_norm[,"Null_M_vs_F_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"Null_M_vs_F_diff"] = log2(fpkm_norm[,"Null_M"] + stab_var) - log2(fpkm_norm[,"Null_F"] + stab_var)

set1=c("Wt_M_1","Wt_M_2","Wt_M_3"); set2=c("Null_M_1","Null_M_2","Null_M_3")
fpkm_norm[,"M_Wt_vs_Null_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"M_Wt_vs_Null_diff"] = log2(fpkm_norm[,"Wt_M"] + stab_var) - log2(fpkm_norm[,"Null_M"] + stab_var)

set1=c("Wt_F_1","Wt_F_2","Wt_F_3"); set2=c("Null_F_1","Null_F_2","Null_F_3")
fpkm_norm[,"F_Wt_vs_Null_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"F_Wt_vs_Null_diff"] = log2(fpkm_norm[,"Wt_F"] + stab_var) - log2(fpkm_norm[,"Null_F"] + stab_var)

set1=c("Wt_M_1","Wt_M_2","Wt_M_3","Null_M_1","Null_M_2","Null_M_3"); set2=c("Wt_F_1","Wt_F_2","Wt_F_3","Null_F_1","Null_F_2","Null_F_3")
fpkm_norm[,"M_vs_F_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"M_vs_F_diff"] = log2(fpkm_norm[,"M"] + stab_var) - log2(fpkm_norm[,"F"] + stab_var)

set1=c("Wt_M_1","Wt_M_2","Wt_M_3","Wt_F_1","Wt_F_2","Wt_F_3"); set2=c("Null_M_1","Null_M_2","Null_M_3","Null_F_1","Null_F_2","Null_F_3")
fpkm_norm[,"Wt_vs_Null_pval"] = apply(fpkm_norm, 1, ttest)
fpkm_norm[,"Wt_vs_Null_diff"] = log2(fpkm_norm[,"Wt"] + stab_var) - log2(fpkm_norm[,"Null"] + stab_var)

#Get final gene lists
data2 = data.frame(data[,c("tracking_id","gene_names")], fpkm_norm)

M_vs_F_i=which(fpkm_norm[,"M_vs_F_pval"] < 0.05 & abs(fpkm_norm[,"M_vs_F_diff"]) > 0.5)
M_vs_F_list=data2[M_vs_F_i,c("gene_names","M_vs_F_diff","M_vs_F_pval")]
write.table(M_vs_F_list, "Male_vs_Female_DE_Genes.txt", sep="\t", row.names=FALSE)

Wt_vs_Null_i=which(fpkm_norm[,"Wt_vs_Null_pval"] < 0.05 & abs(fpkm_norm[,"Wt_vs_Null_diff"]) > 0.5)
Wt_vs_Null_list=data2[Wt_vs_Null_i,c("gene_names","Wt_vs_Null_diff","Wt_vs_Null_pval")]
write.table(Wt_vs_Null_list, "WildType_vs_Null_DE_Genes.txt", sep="\t", row.names=FALSE)

write.table(data2, "All_Data.txt", sep="\t", row.names=FALSE)

####################################################
#Plot all within lot comparisons.  Include 0s.
plotComp2=function(indata, s1, s1_name, s2, s2_name, main_name, comp_name, size){
	#Get min non-zero value from the two vectors
	z=indata
	z[indata==0]=NA
	min1=min(z[,s1], na.rm=TRUE)	
	min2=min(z[,s2], na.rm=TRUE)	
	min_nonzero=min(c(min1,min2))
	min_nonzero=0.0001
	legend_lwd=2
	sigcol="magenta"
	legend_cols=c("black", "orange", "red", sigcol)
	legend_lty=2
	plot(x=log2(indata[,s1]+min_nonzero), y=log2(indata[,s2]+min_nonzero), main=main_name, 
		 xlab=s1_name, ylab=s2_name, col="grey20")
	i=which(fpkm_norm[,comp_name] < 0.05)
	n=length(i)
	points(x=log2(indata[i,s1]+min_nonzero), y=log2(indata[i,s2]+min_nonzero), col=sigcol, pch=20, cex=size)
	abline(0,1, col="black", lwd=2)
	abline(1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
	abline(-1,1, col="orange", lty=legend_lty, lwd=legend_lwd)
	abline(2,1, col="red", lty=legend_lty, lwd=legend_lwd)
	abline(-2,1, col="red", lty=legend_lty, lwd=legend_lwd)
	legend_text = c("No change","Fold change > 2", "Fold change > 4", paste("Significant (n = ", n, ")", sep=""))
	legend("topleft", legend_text, lwd=legend_lwd, lty=c(1,legend_lty,legend_lty,NA), pch=c(NA,NA,NA,20), col=legend_cols)
}

plotComp2(fpkm_norm, "Wt_M", "FPKM (WT, Male, Mean)", "Wt_F", "FPKM (Wt, Female, Mean)", "Correlation between FPKMs, WT Male vs Female, Means", "Wt_M_vs_F_pval", 1)
plotComp2(fpkm_norm, "Null_M", "FPKM (Null, Male, Mean)", "Null_F", "FPKM (Null, Female, Mean)", "Correlation between FPKMs, NULL Male vs Female, Means", "Null_M_vs_F_pval", 1)
plotComp2(fpkm_norm, "Wt_M", "FPKM (WT, Male, Mean)", "Null_M", "FPKM (Null, Male, Mean)", "Correlation between FPKMs, Male WT vs NULL, Means", "M_Wt_vs_Null_pval", 0.75)
plotComp2(fpkm_norm, "Wt_F", "FPKM (WT, Female, Mean)", "Null_F", "FPKM (Null, Female, Mean)", "Correlation between FPKMs, Female WT vs NULL, Means", "F_Wt_vs_Null_pval", 0.5)

plotComp2(fpkm_norm, "M", "FPKM (Male, Mean)", "F", "FPKM (Female, Mean)", "Correlation between FPKMs, Male vs Female, Means", "M_vs_F_pval", 0.5)
plotComp2(fpkm_norm, "Wt", "FPKM (Wt, Mean)", "Null", "FPKM (Null, Mean)", "Correlation between FPKMs, Wt vs Null, Means", "Wt_vs_Null_pval", 0.25)



#Try the MDS plots on a lot-by-lot basis
min_cov=1
mmfilt=function(minv=min_variance){
	function(x){
		cov=sd(x, na.rm=TRUE)/mean(x,na.rm=TRUE)
		(cov > min_cov)
	}	
}
mmfun = mmfilt()
ffun = filterfun(mmfun)
i = genefilter(fpkm, ffun)
sum(i)

getMDS=function(x, i){
	r=cor(x[i,], use="pairwise.complete.obs", method="spearman")
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	return(mds)
}
lot1=c(1,4,7,10)
lot2=c(2,5,8,11)
lot3=c(3,6,9,12)
mds1=getMDS(fpkm[,lot1],i)
mds2=getMDS(fpkm[,lot2],i)
mds3=getMDS(fpkm[,lot3],i)
plotMDS=function(mds, lot, label){
	zz=0.2
	min_x=min(mds$points[,1])-abs((min(mds$points[,1])*zz))
	max_x=max(mds$points[,1])+abs((max(mds$points[,1])*zz))
	min_y=	min(mds$points[,2])-abs((min(mds$points[,2])*zz))
	max_y=max(mds$points[,2])+abs((max(mds$points[,2])*zz))
	main=paste("(cov > ", min_cov, "; n genes = ", sum(i), ")", sep="")
	plot(mds$points, type="n", xlab="", ylab="", main=paste(label, ": Classes", main), xlim=c(min_x,max_x), ylim=c(min_y,max_y)); text(mds$points[,1], mds$points[,2], full_class[lot], col=class_cols[lot])
}
par(mfrow=c(2,2))
plotMDS(mds1, lot1, "Lot1")
plotMDS(mds2, lot2, "Lot2")
plotMDS(mds3, lot3, "Lot3")

#Sanity check of NF1 expression (I determined which data line corresponds to NF1 by grepping the following line and then matching the coords/exp values)
#/gscmnt/sata132/techd/solexa/jwalker/RNAseq/Mouse/Josh_Rubin/Josh_Rubin_ab_initio/cuffdiff_wGene_multiExon/genes.fpkm_tracking
nf1_id="XLOC_017620"
nf1_i=which(data[,"tracking_id"]==nf1_id)

n=list(c("Wt","Null"), c("Male1","Male2","Male3","Female1","Female2","Female3"))
nf1_fpkm_mat=matrix(t(fpkm[nf1_i,c(1,7,2,8,3,9,4,10,5,11,6,12)]), ncol=6, dimnames=n)
barplot(nf1_fpkm_mat, beside=TRUE, col=c("blue","red"), ylab="NF1 FPKM", xlab="Sample pair (red=WT, blue=Null)", main="NF1 expression in WT vs. Null")

nf1_fpkm_mat=matrix(t(fpkm_norm[nf1_i,c(1,7,2,8,3,9,4,10,5,11,6,12)]), ncol=6, dimnames=n)
barplot(nf1_fpkm_mat, beside=TRUE, col=c("blue","red"), ylab="NF1 FPKM (Normalized)", xlab="Sample pair (red=WT, blue=Null)", main="NF1 expression in WT vs. Null")

#Expression of a specific gene
plotGeneExp1 = function(gene_name){
	comp=c(1,4,2,5,3,6,7,10,8,11,9,12)
	gene_i=which(data[,"gene_names"]==gene_name)
	n=list(c("Male","Female"), c("Wt1","Wt2","Wt3","Null1","Null2","Null3"))
	gene_fpkm_mat=matrix(t(fpkm[gene_i,comp]), ncol=6, dimnames=n)
	par(mfrow=c(2,1))
	ylab1=paste(gene_name, " FPKM", sep="")
	ylab2=paste(gene_name, " FPKM (Normalized)", sep="")
	title=paste(gene_name, " expression in male vs. female")
	barplot(gene_fpkm_mat, beside=TRUE, col=c("blue","pink"), ylab=ylab1, xlab="Sample pair (blue=Male, pink=Female)", main=title)
	gene_fpkm_mat=matrix(t(fpkm_norm[gene_i,comp]), ncol=6, dimnames=n)
	barplot(gene_fpkm_mat, beside=TRUE, col=c("blue","pink"), ylab=ylab2, xlab="Sample pair (blue=Male, pink=Female)", main=title)	
}
plotGeneExp2 = function(gene_name){
	comp=c(1,7,2,8,3,9,4,10,5,11,6,12)
	gene_i=which(data[,"gene_names"]==gene_name)
	n=list(c("Wt","Null"), c("Wt1","Wt2","Wt3","Null1","Null2","Null3"))
	gene_fpkm_mat=matrix(t(fpkm[gene_i,comp]), ncol=6, dimnames=n)
	par(mfrow=c(2,1))
	ylab1=paste(gene_name, " FPKM", sep="")
	ylab2=paste(gene_name, " FPKM (Normalized)", sep="")
	title=paste(gene_name, " expression in WT vs. Null")
	barplot(gene_fpkm_mat, beside=TRUE, col=c("blue","pink"), ylab=ylab1, xlab="Sample pair (blue=Male, pink=Female)", main=title)
	gene_fpkm_mat=matrix(t(fpkm_norm[gene_i,comp]), ncol=6, dimnames=n)
	barplot(gene_fpkm_mat, beside=TRUE, col=c("blue","pink"), ylab=ylab2, xlab="Sample pair (blue=Male, pink=Female)", main=title)	
}


#Sanity tests from the literature
#Male vs. Female single gene comparisons

#X chromosome
plotGeneExp1("Utx")

#Y chromosome
plotGeneExp1("Uty")
plotGeneExp1("AK020213,Ddx3y,dby")
plotGeneExp1("Jarid1d,Smcy")
plotGeneExp1("Eif2s3y")


#Wt. vs. Null single gene comparisons
plotGeneExp2("NF1GRP,Nf1,NF1")























