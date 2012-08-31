#Create an MDS plot for the BRAF inhibitor resistant cell lines using their VAFs
data_dir = "/Users/malachig/Dropbox/Documents/Analysis_development/braf"
setwd(data_dir)
dir()

vafs = read.table("BRAF_WGS_SNV_VAFs_Matrix.tsv", header=TRUE, as.is=c(1))
muts = read.table("BRAF_WGS_SNV_MutationStatus_Matrix.tsv", header=TRUE, as.is=c(1))
sample_names = names(vafs)[2:length(vafs)]
sample_names_abr = c("F-R12", "F-R6", "F-R8", "V-R10", "V-R19", "V-R22", "V-R23", "V-1R")
sample_cols = c("blue","blue","blue","dark green","dark green","dark green","dark green","dark green")
sample_pch = c(16,16,16,17,17,17,17,17)

#Define function to create MDS plots
plotMDS = function(filename, data, sample_names, sample_names_abr, cor_method, title, xlimit, ylimit, sample_pch, sample_cols, points, text){
	pdf(filename)
	r=cor(data[,sample_names], use="pairwise.complete.obs", method=cor_method)
	r[is.na(r)]=0
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	xlabel = paste("X MDS distance (", cor_method, ")", sep="")
	ylabel = paste("Y MDS distance (", cor_method, ")", sep="")	
	plot(mds$points, type="n", main=title, xlim=xlimit, ylim=ylimit, xlab=xlabel, ylab=ylabel)
	if(points){
		points(mds$points[,1], mds$points[,2], pch=sample_pch, col=sample_cols)
	}
	if(text){
		text(mds$points[,1], mds$points[,2], sample_names_abr, col=sample_cols)
	}
	dev.off()
}

#Visualize the correlation values
cor(vafs[,sample_names], , use="pairwise.complete.obs", method="pearson")

#Points, pearson, VAFs
plotMDS("MDS_Points_Pearson_VAFs.pdf", vafs, sample_names, sample_names_abr, "pearson", "MDS distance plot using VAFs", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 1, 0)

#Text, pearson, VAFs
plotMDS("MDS_Text_Pearson_VAFs.pdf", vafs, sample_names, sample_names_abr, "pearson", "MDS distance plot using VAFs", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 0, 1)

#Points, spearman, VAFs
plotMDS("MDS_Points_Spearman_VAFs.pdf", vafs, sample_names, sample_names_abr, "spearman", "MDS distance plot using VAFs", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 1, 0)

#Text, spearman, VAFs
plotMDS("MDS_Text_Spearman_VAFs.pdf", vafs, sample_names, sample_names_abr, "spearman", "MDS distance plot using VAFs", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 0, 1)

#Points, pearson, MUTs
plotMDS("MDS_Points_Pearson_MUTs.pdf", muts, sample_names, sample_names_abr, "pearson", "MDS distance plot using mutation status", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 1, 0)

#Text, pearson, MUTs
plotMDS("MDS_Text_Pearson_MUTs.pdf", muts, sample_names, sample_names_abr, "pearson", "MDS distance plot using mutation status", c(-0.7,0.65), c(-0.7,0.6), sample_pch, sample_cols, 0, 1)

#Points, spearman, MUTs
plotMDS("MDS_Points_Spearman_MUTs.pdf", muts, sample_names, sample_names_abr, "spearman", "MDS distance plot using mutation status", c(-0.7,0.6), c(-0.7,0.6), sample_pch, sample_cols, 1, 0)

#Text, spearman, MUTs
plotMDS("MDS_Text_Spearman_MUTs.pdf", muts, sample_names, sample_names_abr, "spearman", "MDS distance plot using mutation status", c(-0.7,0.65), c(-0.7,0.6), sample_pch, sample_cols, 0, 1)







