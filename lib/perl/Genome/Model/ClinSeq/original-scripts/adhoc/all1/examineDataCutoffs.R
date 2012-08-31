#Pre-filtered data
datadir="/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/result/pre_filtered/"
outdir = "/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/figures/pre_filtered"

#Post-filtered data
datadir="/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/result/post_filtered/"
outdir = "/gscmnt/sata132/techd/mgriffit/hgs/custom/all01/snp_array_design/figures/post_filtered/"

setwd(datadir)
data=read.table("final_target_positions.tsv", header=TRUE, sep="\t", as.is=c(1,4:11,20))
dir()
setwd(outdir)
dir()


#READ COVERAGE - Plot coverage distribution from the Tier1 SNVs only and use this distribution to calculate a cutoff for all SNP positions
#TIER 1 SNVS
min_normal_read_coverage = 20
i=which(data[,"class"] == "Tier1_SNV_Novel")
pdf("Tier1_SNV_ReadCoverage_Boxplot.pdf")
b=boxplot(data[i,"normal_cov"], col="red", main="Box plot of read coverage for all Tier1 SNV positions", ylab="Tier1 SNV read coverage values (WGS + Exome reads)", xlab="Tier1 SNVs")
max_normal_read_coverage=b$stats[5,1]
abline(h=max_normal_read_coverage, lty=2, lwd=2)
abline(h=min_normal_read_coverage, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()
b$stats[5,1] #Get the upper cutoff for outliers from this boxplot

#ALL SNVS
i=which(data[,"class"] == "Tier1_SNV_Known" | data[,"class"] == "Tier1_SNV_Novel" | data[,"class"] == "Tier2_SNV_Known" | data[,"class"] == "Tier2_SNV_Novel" | data[,"class"] == "Tier3_SNV_Known" | data[,"class"] == "Tier3_SNV_Novel")
pdf("ReadCoverage_Tier1-3_SNV_Boxplot.pdf")
b=boxplot(data[i,"normal_cov"], col="red", main="Box plot of read coverage for Tier1-3 SNV positions", ylab="Tier1-3 SNV read coverage values (WGS + Exome reads)", xlab="Tier1-3 SNVs")
abline(h=max_normal_read_coverage, lty=2, lwd=2)
abline(h=min_normal_read_coverage, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#ALL POSITIONS (SNVS + SNPS)
pdf("ReadCoverage_AllPositions_Boxplot.pdf")
b=boxplot(data[,"normal_cov"], col="red", main="Box plot of read coverage for all SNV+SNP positions", ylab="Normal read coverage values (WGS + Exome reads)", xlab="All SNVs and SNPs")
abline(h=max_normal_read_coverage, lty=2, lwd=2)
abline(h=min_normal_read_coverage, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#DELETION SNP POSITIONS ONLY
i = which(data[,"class"] == "Deletion_Het_SNP")
pdf("ReadCoverage_DeletionPositions_Boxplot.pdf")
b=boxplot(data[i,"normal_cov"], col="red", main="Box plot of read coverage for deletion positions", ylab="Normal read coverage values (WGS + Exome reads)", xlab="Deletion SNPs")
abline(h=max_normal_read_coverage, lty=2, lwd=2)
abline(h=min_normal_read_coverage, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#CONTROL SNP POSITIONS ONLY
i = which(data[,"class"] == "Control_Het_SNP")
pdf("ReadCoverage_ControlPositions_Boxplot.pdf")
b=boxplot(data[i,"normal_cov"], col="red", main="Box plot of read coverage for all control positions", ylab="Normal read coverage values (WGS + Exome reads)", xlab="Control SNPs")
abline(h=max_normal_read_coverage, lty=2, lwd=2)
abline(h=min_normal_read_coverage, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#Read coverage by position type
pdf("ReadCoverage_ByPositionType_Boxplot.pdf")
op=par()
par(mar=c(10, 4, 4, 2)+0.1)
boxplot(log2(data[,"normal_cov"]+1)~data[,"class"], main="Box plot of read coverage for all SNV+SNP positions", ylab="Log2 SNV+SNP read coverage values (WGS + Exome reads)", xlab="", las=2, col=rainbow(length(unique(data[,"class"]))), varwidth=TRUE)
par(op)
dev.off()


#NORMAL VAFS - Plot distribution of 'normal_vaf'
#ALL POSITIONS
pdf("Normal_VAFs_AllPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[,"normal_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Normal* Sample", ylab="Density", main="Distribution of Normal VAFs for *all selected* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#ALL SNVS
i=which(data[,"class"] == "Tier1_SNV_Known" | data[,"class"] == "Tier1_SNV_Novel" | data[,"class"] == "Tier2_SNV_Known" | data[,"class"] == "Tier2_SNV_Novel" | data[,"class"] == "Tier3_SNV_Known" | data[,"class"] == "Tier3_SNV_Novel")
pdf("Normal_VAFs_Tier1-3_SnvPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"normal_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Normal* Sample", ylab="Density", main="Distribution of Normal VAFs for *all SNV* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#DELETION SNP POSITIONS ONLY
i = which(data[,"class"] == "Deletion_Het_SNP")
pdf("Normal_VAFs_DeletionPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"normal_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Normal* Sample", ylab="Density", main="Distribution of Normal VAFs for *deletion* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#CONTROL SNP POSITIONS ONLY
i = which(data[,"class"] == "Control_Het_SNP")
pdf("Normal_VAFs_ControlPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"normal_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Normal* Sample", ylab="Density", main="Distribution of Normal VAFs for *control* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()


#TUMOR VAFS - Plot distribution of 'tumor_vaf'
#ALL POSITIONS
pdf("Tumor_VAFs_AllPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[,"tumor_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Tumor* Sample", ylab="Density", main="Distribution of Tumor VAFs for *all selected* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#ALL SNVS
i=which(data[,"class"] == "Tier1_SNV_Known" | data[,"class"] == "Tier1_SNV_Novel" | data[,"class"] == "Tier2_SNV_Known" | data[,"class"] == "Tier2_SNV_Novel" | data[,"class"] == "Tier3_SNV_Known" | data[,"class"] == "Tier3_SNV_Novel")
pdf("Tumor_VAFs_Tier1-3_SnvPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"tumor_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Tumor* Sample", ylab="Density", main="Distribution of Tumor VAFs for *all SNV* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#DELETION SNP POSITIONS ONLY
i = which(data[,"class"] == "Deletion_Het_SNP")
pdf("Tumor_VAFs_DeletionPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"tumor_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Tumor* Sample", ylab="Density", main="Distribution of Tumor VAFs for *deletion* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#CONTROL SNP POSITIONS ONLY
i = which(data[,"class"] == "Control_Het_SNP")
pdf("Tumor_VAFs_ControlPositions_Hist.pdf")
het_vaf_diff_var = 10
hist(data[i,"tumor_vaf"], breaks=50, col="purple2", xlab="Variant Allele Frequency (VAF) in *Tumor* Sample", ylab="Density", main="Distribution of Tumor VAFs for *control* positions")
abline(v=50-het_vaf_diff_var, lty=2, lwd=2)
abline(v=50+het_vaf_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()


#CNV WINDOW DIFFERENCE VALUES - Plot distribution of 'cnv_window_diff' for deletion regions only
#DELETION SNP POSITIONS ONLY
cnv_deletion_window_diff_var = 0 - 0.25
i = which(data[,"class"] == "Deletion_Het_SNP")
pdf("CNV_Diffs_DeletionPositionsOnly_Hist.pdf")
hist(data[i,"cnv_window_diff"], breaks=50, xlim=c(-2.5,1), col="dark blue", xlab="CNV difference between tumor and normal", ylab="Density", main="CNV diffs for SNP positions within deletion regions")
abline(v=cnv_deletion_window_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#CONTROL SNP POSITIONS ONLY
i = which(data[,"class"] == "Control_Het_SNP")
cnv_control_window_diff_var = 0.10
pdf("CNV_Diffs_ControlPositionsOnly_Hist.pdf")
hist(data[i,"cnv_window_diff"], breaks=250, xlim=c(-2.5,1), col="blue", xlab="CNV difference between tumor and normal", ylab="Density", main="CNV diffs for SNP positions within control regions")
abline(v=0-cnv_control_window_diff_var, lty=2, lwd=2)
abline(v=0+cnv_control_window_diff_var, lty=2, lwd=2)
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()

#TIER 1-3 SNVS ONLY
i=which(data[,"class"] == "Tier1_SNV_Known" | data[,"class"] == "Tier1_SNV_Novel" | data[,"class"] == "Tier2_SNV_Known" | data[,"class"] == "Tier2_SNV_Novel" | data[,"class"] == "Tier3_SNV_Known" | data[,"class"] == "Tier3_SNV_Novel")
pdf("CNV_Diffs_SnvPositionsOnly_Hist.pdf")
hist(data[i,"cnv_window_diff"], breaks=250, xlim=c(-2.5,1), col="blue", xlab="CNV difference between tumor and normal", ylab="Density", main="CNV diffs for all SNV positions")
legend("topleft", "Cutoff", lty=2, lwd=2, bg="white")
dev.off()


#CLONALITY PLOTS - TUMOR VAFS VS. READ COVERAGE
#TIER1 SNVS ONLY
#Look at the variant allele frequencies vs. read coverage for the Tier1 SNVs to allow comparison back to the clonality analysis plot
#tumor_cov vs. tumor_vaf
i=which(data[,"class"] == "Tier1_SNV_Novel")
pchs = rep(16, length(data[i,"chr"]))
pchs[which(data[i,"chr"]=="X")] = 17
minor_clone_vaf=20.2
major_clone_vaf=49.7
pdf("Clonality_Tier1_SNVs_Only.pdf")
plot(x=data[i,"tumor_vaf"], y=data[i,"tumor_cov"], pch=pchs, col="light green", xlab="Tumor VAF", ylab="Tumor Read Coverage", main="ALL1 - Tumor Clonality - Tier1 SNVs only")
abline(v=minor_clone_vaf, lty=2, lwd=2)
abline(v=major_clone_vaf, lty=3, lwd=2)
legend("topleft", c("Minor clone VAF", "Major clone VAF"), lty=c(2,3), lwd=2, bg="white")
dev.off()

#TIER 1-3 SNVS ONLY
i=which(data[,"class"] == "Tier1_SNV_Known" | data[,"class"] == "Tier1_SNV_Novel" | data[,"class"] == "Tier2_SNV_Known" | data[,"class"] == "Tier2_SNV_Novel" | data[,"class"] == "Tier3_SNV_Known" | data[,"class"] == "Tier3_SNV_Novel")
pchs = rep(16, length(data[i,"chr"]))
pchs[which(data[i,"chr"]=="X")] = 17
pcols = rep("light green", length(data[i,"chr"]))
pcols[which(data[i,"chr"]=="X")] = "orange"
minor_clone_vaf=20.2
major_clone_vaf=49.7
pdf("Clonality_Tier1-3_SNVs_Only.pdf")
plot(x=data[i,"tumor_vaf"], y=data[i,"tumor_cov"], pch=pchs, col=pcols, xlab="Tumor VAF", ylab="Tumor Read Coverage", main="ALL1 - Tumor Clonality - Tier1 SNVs only", cex=0.35)
abline(v=minor_clone_vaf, lty=2, lwd=2)
abline(v=major_clone_vaf, lty=3, lwd=2)
legend("topleft", c("Minor clone VAF", "Major clone VAF"), lty=c(2,3), lwd=2, bg="white")
dev.off()


#Pie chart of position classes
pdf("ClassTypes_Pie.pdf")
pie(table(data[,"class"]), col=rainbow(length(unique(data[,"class"]))))
dev.off()

#Pie chart of variant types
pdf("VariantTypes_Pie.pdf")
pie(table(data[,"var_type"]), col=rainbow(length(unique(data[,"class"]))))
dev.off()


