
#datadir="/Users/mgriffit/Desktop/HG1/RNA-seq"
datadir="/gscmnt/sata132/techd/mgriffit/hg1/rna_seq/hg1_vs_brcs/"
setwd(datadir)
dir()

library(preprocessCore)
header=c("tracking_id","class_code","nearest_ref_id", "gene_id", "gene_short_name", "tss_id", "locus", "length", "coverage", "status", "FPKM", "FPKM_conf_lo", "FPKM_conf_hi")
data_brc4 = read.table("BRC4.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_brc4) = header
data_brc5 = read.table("BRC5.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_brc5) = header
data_brc18 = read.table("BRC18.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_brc18) = header
data_brc31 = read.table("BRC31.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_brc31) = header
data_hg1 = read.table("HG1.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_hg1) = header
data_hg1_lib1 = read.table("HG1_cDNA1.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_hg1_lib1) = header
data_hg1_lib2 = read.table("HG1_cDNA2.genes.fpkm_tracking.sort", header=FALSE, as.is=1:10)
names(data_hg1_lib2) = header

#Just a bar plot of ERBB2 by itself using raw data
gene_name="ERBB2"
gene_name="ESR1"
gene_name="PGR"
gene_name="HDAC2"
i=which(data_brc4[,"gene_short_name"]==gene_name)


fpkm=data.frame(data_brc4[,"FPKM"], data_brc5[,"FPKM"], data_brc18[,"FPKM"], data_brc31[,"FPKM"], data_hg1[,"FPKM"])
names(fpkm) = c("BRC4", "BRC5", "BRC18", "BRC31", "HG1")
fpkm_low=data.frame(data_brc4[,"FPKM_conf_lo"], data_brc5[,"FPKM_conf_lo"], data_brc18[,"FPKM_conf_lo"], data_brc31[,"FPKM_conf_lo"], data_hg1[,"FPKM_conf_lo"])
fpkm_hi=data.frame(data_brc4[,"FPKM_conf_hi"], data_brc5[,"FPKM_conf_hi"], data_brc18[,"FPKM_conf_hi"], data_brc31[,"FPKM_conf_hi"], data_hg1[,"FPKM_conf_hi"])

main_lab = paste(gene_name, " - RNA-seq expression level (Cufflinks)", sep="")
max_y=max(as.numeric(fpkm[i,]))
upper_y = max_y+(max_y*0.2)
bp=barplot(as.numeric(fpkm[i,]), col=c("honeydew","honeydew","honeydew","honeydew","turquoise3"), ylim=c(0,upper_y), names=names(fpkm), xlab="Library", ylab="FPKM", main=main_lab)
arrows(x0=bp, y0=as.numeric(fpkm_low[i,]), x1=bp, y1=as.numeric(fpkm_hi[i,]), length=0.15, angle=90, code=3, lwd=2)

#Now normalize the data and plot ERBB2 on the distribution of gene expression for each library
names=names(fpkm)
x = as.matrix(fpkm)
fpkm_norm = normalize.quantiles(x)
fpkm_norm = as.data.frame(fpkm_norm)
dimnames(fpkm_norm)[[2]] = names
fpkm=fpkm_norm

#Create a box plot of all distributions
z=0.00000001
fpkm2=fpkm
fpkm2[fpkm==0]=NA
fpkm=log2(fpkm2+z)
 
names(fpkm) = c("BRC4", "BRC5", "BRC18", "BRC31", "HG1")
x=list(fpkm[,"BRC4"], fpkm[,"BRC5"], fpkm[,"BRC18"], fpkm[,"BRC31"], fpkm[,"HG1"])
names(x) = c("BRC4", "BRC5", "BRC18", "BRC31", "HG1")
main_lab = paste(gene_name, " - Distributions of all gene FPKMs (zeros removed)", sep="")
bp=boxplot(x, col=c("honeydew","honeydew","honeydew","honeydew","turquoise3"), xlab="Library", ylab="Log2 FPKM", main=main_lab)

#Mark the expression level of ERBB2 on each
points(x=1:5, y=fpkm[i,], col="red", pch=16)
legend("topleft", gene_name, col="red", pch=16, bg="white")


#Evaluate the HG1 RNA-seq 'replicates' 
data_hg1_comp = data.frame(data_hg1_lib1[,"FPKM"], data_hg1_lib2[,"FPKM"])
names(data_hg1_comp) = c("FPKM_Lib1", "FPKM_Lib2")
x = as.matrix(data_hg1_comp)
x_norm = normalize.quantiles(x)
x_norm = as.data.frame(x_norm)
dimnames(x_norm)[[2]] = c("FPKM_Lib1", "FPKM_Lib2")
data_hg1_comp=x_norm

x=log2(data_hg1_comp[,"FPKM_Lib1"]+1)
y=log2(data_hg1_comp[,"FPKM_Lib2"]+1)

#Plot all the data
main_lab="Expression correlation (library 1 vs 2)"
x_lab="Library 1 - log2(FPKM+1)"
y_lab="Library 2 - log2(FPKM+1)"
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jpeg("Lib1_vs_Lib2.jpg", quality=100, type="quartz", width=600, height=600)
smoothScatter(x=x, y=y, main=main_lab, xlab=x_lab, ylab=y_lab, colramp=colors, nbin=275,
                col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()

#Apply an expression level cutoff
data_hg1_comp[,"FPKM_SUM"] = data_hg1_comp[,"FPKM_Lib1"] + data_hg1_comp[,"FPKM_Lib2"]
i=which(data_hg1_comp[,"FPKM_SUM"] > 2)
jpeg("Lib1_vs_Lib2_Cutoff.jpg", quality=100, type="quartz", width=600, height=600)
smoothScatter(x=x[i], y=y[i], main=main_lab, xlab=x_lab, ylab=y_lab, colramp=colors, nbin=275,
                col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()


#Calculate a simple fold change for all BRCs vs. HG1
fpkm=data.frame(data_brc4[,"FPKM"], data_brc5[,"FPKM"], data_brc18[,"FPKM"], data_brc31[,"FPKM"], data_hg1[,"FPKM"])
sample_names = c("BRC4", "BRC5", "BRC18", "BRC31", "HG1", "SYMBOL") 
x = as.matrix(fpkm)
x_norm = normalize.quantiles(x)
x_norm = as.data.frame(x_norm)
fpkm=x_norm
fpkm[,"SYMBOL"] = data_brc4[,"gene_short_name"]
names(fpkm) = sample_names

fpkm[,"BRC_MEAN"] = apply(fpkm[,c("BRC4", "BRC5", "BRC18", "BRC31")], 1, mean)
fpkm[,"FPKM_DIFF"] = fpkm[,"HG1"] - fpkm[,"BRC_MEAN"] 
o = order(abs(fpkm[,"FPKM_DIFF"]), decreasing=TRUE)
write.table(fpkm[o, ],"HG1_vs_BRCs_FPKM_Diff.tsv", sep="\t", row.names=FALSE, quote=FALSE)

x=log2(fpkm[,"BRC_MEAN"]+1)
y=log2(fpkm[,"HG1"]+1)

jpeg("HG1_vs_BRCs_FPKM_DE.jpg", quality=100, type="quartz", width=600, height=600)
smoothScatter(x=x, y=y, main="FPKM Difference (HG1 vs BRCs)", xlab="BRC Mean FPKM", ylab="HG1 FPKM", colramp=colors, nbin=275,
                col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
abline(0, 1, lty=2)
dev.off()

fpkm[o[1:100],]
