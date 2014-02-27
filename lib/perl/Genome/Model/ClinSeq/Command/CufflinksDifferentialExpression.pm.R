#!/usr/bin/env Rscript
#Written by Malachi Griffith

#- produce de file for all genes with fold change and de status (de hq, de lq, no change)
#- determine fold change as log2 difference?
#- hq will be those diffs where the conf intervals do not overlap and if available, both have status of 'OK'
#- produce .hq and .lq differential expression files for both gains, losses, and gains+losses
#- 'de' will be those that exceed fold-change cutoff
#- files:
#- case_vs_control.tsv
#- case_vs_control.hq.de.tsv, case_vs_control.hq.up.tsv, case_vs_control.hq.down.tsv
#- case_vs_control.lq.de.tsv, case_vs_control.lq.up.tsv, case_vs_control.lq.down.tsv
#- Create a plot showing the distribution of FPKM values for both samples
#- Create a plot of FPKM case vs. FPKM control.  Color plot with hq de genes/transcripts

args = (commandArgs(TRUE))
outdir = args[1];        #working directory -> used to load input files and dump results files
infile = args[2];        #gene.de.input.tsv or transcript.de.input.tsv -> input file containing case and control FPKM values
type = args[3];          #'gene' or 'transcript' -> used for naming output files
case_label = args[4];    #e.g. 'tumor'
control_label = args[5]; #e.g. 'normal'

if (length(args) < 2){
  message_text1 = "Required arguments missing: ./CufflinksDifferentialExpression.pm.R /tmp/rnaseq_de gene.de.input.tsv gene tumor normal"
  stop(message_text1)
}
#outdir = "/tmp/rnaseq_de";
#infile = "gene.de.input.tsv";
#type = "gene"
#case_label = "tumor"
#control_label = "normal"

#load libraries
library(preprocessCore)
library(ggplot2)


#Set global parameters that determine genes/transcripts that will be considered DE
variance_stabilization = 0.1
min_diff = 1

setwd(outdir)
data = read.table(file=infile, sep="\t", header=TRUE, as.is=c(1:4), na.strings = c("NA", "na"))

#Create some plots to summarize the data pre normalization

#Scatter plot of case vs. control FPKM values before normalization
main_label = paste(control_label, " vs. ", case_label, " expression (", type, ")", sep="")
x_label = paste(case_label, " - Cufflinks FPKM (log2)", sep="")
y_label = paste(control_label, " - Cufflinks FPKM (log2)", sep="")
png(filename = "case_vs_control_fpkm_scatter_prenorm.png", width = 800, height = 800, bg = "white")
plot(x=log2(data[,"case_fpkm"]+variance_stabilization), y=log2(data[,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="blue", xlab=x_label, ylab=y_label, main=main_label)
abline(a=0, b=1, lwd=2, lty=2, col="black")
abline(a=min_diff, b=1, lwd=1, lty=2, col="black")
abline(a=min_diff*-1, b=1, lwd=1, lty=2, col="black")
dev.off()

#Density plot of case vs. control FPKM values before normalization
x = data[which(data[,"case_fpkm"] > 0), "case_fpkm"]
y = data[which(data[,"control_fpkm"] > 0), "control_fpkm"]
samples = c(rep(case_label, length(x)), rep(control_label, length(y)))
fpkms = data.frame(c(log2(x+variance_stabilization),log2(y+variance_stabilization)), samples)
names(fpkms) = c("FPKM","Sample")
fpkm_count1 = length(x)
fpkm_count2 = length(y)
y_label = paste("Density (n = ", fpkm_count1, " and ", fpkm_count2, " ", type, "s)", sep="")
pdf("case_vs_control_fpkm_density_prenorm.pdf")
if (packageVersion("ggplot2") <= "0.8.9"){
  print ({
    ggplot(fpkms, aes(FPKM, fill = Sample)) + geom_density(alpha = 0.2) + xlab("FPKM (log2)") + ylab(y_label) + opts(title=main_label)
  })
}else{
  print ({
    ggplot(fpkms, aes(FPKM, fill = Sample)) + geom_density(alpha = 0.2) + xlab("FPKM (log2)") + ylab(y_label) + labs(title=main_label)
  })
}
dev.off()


#Perform a normalization of the fpkm values
x1 = data[,c("case_fpkm","control_fpkm")]
x2 = as.matrix(x1)
x3 = normalize.quantiles(x2)
x3 = as.data.frame(x3)
data[,"case_fpkm"] = x3[,1]
data[,"control_fpkm"] = x3[,2]

#Perform a normalization of the fpkm_conf_hi values
x1 = data[,c("case_fpkm_conf_hi","control_fpkm_conf_hi")]
x2 = as.matrix(x1)
x3 = normalize.quantiles(x2)
x3 = as.data.frame(x3)
data[,"case_fpkm_conf_hi"] = x3[,1]
data[,"control_fpkm_conf_hi"] = x3[,2]

#Perform a normalization of the fpkm_conf_lo values
x1 = data[,c("case_fpkm_conf_lo","control_fpkm_conf_lo")]
x2 = as.matrix(x1)
x3 = normalize.quantiles(x2)
x3 = as.data.frame(x3)
data[,"case_fpkm_conf_lo"] = x3[,1]
data[,"control_fpkm_conf_lo"] = x3[,2]

#Scatter plot of case vs. control FPKM values after normalization
main_label = paste(control_label, " vs. ", case_label, " expression (", type, ")", sep="")
x_label = paste(case_label, " - Cufflinks FPKM (log2)", sep="")
y_label = paste(control_label, " - Cufflinks FPKM (log2)", sep="")
png(filename = "case_vs_control_fpkm_scatter_postnorm.png", width = 800, height = 800, bg = "white")
plot(x=log2(data[,"case_fpkm"]+variance_stabilization), y=log2(data[,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="blue", xlab=x_label, ylab=y_label, main=main_label)
abline(a=0, b=1, lwd=2, lty=2, col="black")
abline(a=min_diff, b=1, lwd=1, lty=2, col="black")
abline(a=min_diff*-1, b=1, lwd=1, lty=2, col="black")
dev.off()

#Density plot of case vs. control FPKM values after normalization
x = data[which(data[,"case_fpkm"] > 0), "case_fpkm"]
y = data[which(data[,"control_fpkm"] > 0), "control_fpkm"]
samples = c(rep(case_label, length(x)), rep(control_label, length(y)))
fpkms = data.frame(c(log2(x+variance_stabilization),log2(y+variance_stabilization)), samples)
names(fpkms) = c("FPKM","Sample")
fpkm_count1 = length(x)
fpkm_count2 = length(y)
y_label = paste("Density (n = ", fpkm_count1, " and ", fpkm_count2, " ", type, "s)", sep="")
pdf("case_vs_control_fpkm_density_postnorm.pdf")
if (packageVersion("ggplot2") <= "0.8.9"){
  print ({
    ggplot(fpkms, aes(FPKM, fill = Sample)) + geom_density(alpha = 0.2) + xlab("FPKM (log2)") + ylab(y_label) + opts(title=main_label)
  })
}else{
  print ({
    ggplot(fpkms, aes(FPKM, fill = Sample)) + geom_density(alpha = 0.2) + xlab("FPKM (log2)") + ylab(y_label) + labs(title=main_label)
  })
}
dev.off()


#Identify hq and lq DE events
data[,"case_vs_control_log2_de"] = log2((data[,"case_fpkm"])+variance_stabilization) - log2((data[,"control_fpkm"])+variance_stabilization);
data[,"case_vs_control_de_lq"] = 0
data[,"case_vs_control_de_hq"] = 0

de_lq_i = which(abs(data[,"case_vs_control_log2_de"]) > min_diff)
data[de_lq_i,"case_vs_control_de_lq"] = 1


#de_hq_i = which((abs(data[,"case_vs_control_log2_de"]) > min_diff) & (data[,"case_fpkm_status"] == "OK") & (data[,"control_fpkm_status"] == "OK") & ((data[,"case_fpkm_conf_hi"] <= data[,"control_fpkm_conf_lo"]) | (data[,"case_fpkm_conf_lo"] >= data[,"control_fpkm_conf_hi"])))
de_hq_i = which((abs(data[,"case_vs_control_log2_de"]) > min_diff) & (data[,"case_fpkm_status"] == "OK") & (data[,"control_fpkm_status"] == "OK") & ((data[,"case_fpkm_conf_hi"] < data[,"control_fpkm_conf_lo"]) | (data[,"case_fpkm_conf_lo"] > data[,"control_fpkm_conf_hi"])))

data[de_hq_i,"case_vs_control_de_hq"] = 1


#Create some plots to summarize the final differential expression results

#Plot histogram of log2 DE values
de_data = data[,c("case_fpkm","control_fpkm","case_vs_control_log2_de")]
title_text = paste("Distribution of DE values (", type, ")", sep="")
x_label = paste("Log2 Differential Expression", sep="")
y_label = paste("Density (n = ", length(de_data[,1]), " ", type,"s)", sep="")
pdf("case_vs_control_de_hist.pdf")
if (packageVersion("ggplot2") <= "0.8.9"){
  print({
    m <- ggplot(de_data, aes(x=case_vs_control_log2_de)); m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density() + 
         opts(title=title_text) + xlab(x_label) + ylab(y_label)
  })
}else{
  print({
    m <- ggplot(de_data, aes(x=case_vs_control_log2_de)); m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density() + 
         labs(title=title_text) + xlab(x_label) + ylab(y_label)
  })
}
dev.off();

#Remove 0 vs. 0 comparisons and replot histogram
i = which(de_data[,"case_fpkm"] > 0 & de_data[,"control_fpkm"] > 0)
de_data = de_data[i,]
title_text = paste("Distribution of DE values (", type, ")", sep="")
x_label = paste("Log2 Differential Expression (0 vs 0 data points removed)", sep="")
y_label = paste("Density (n = ", length(de_data[,1]), " ", type,"s)", sep="")
pdf("case_vs_control_de_filtered_hist.pdf")
if (packageVersion("ggplot2") <= "0.8.9"){
  print({
    m <- ggplot(de_data, aes(x=case_vs_control_log2_de)); m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density() + 
         opts(title=title_text) + xlab(x_label) + ylab(y_label)
  })
}else{
  print({
    m <- ggplot(de_data, aes(x=case_vs_control_log2_de)); m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density() + 
         labs(title=title_text) + xlab(x_label) + ylab(y_label)
  })
}
dev.off();

#Scatter plot of case vs. control FPKM values after normalization with hq and lq features indicated
main_label = paste(case_label, " vs. ", control_label, " expression (", type, ")", sep="")
x_label = paste(case_label, " - Cufflinks FPKM (log2)", sep="")
y_label = paste(control_label, " - Cufflinks FPKM (log2)", sep="")
png(filename = "case_vs_control_fpkm_scatter_postnorm_de.png", width = 800, height = 800, bg = "white")
plot(x=log2(data[,"case_fpkm"]+variance_stabilization), y=log2(data[,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="blue", xlab=x_label, ylab=y_label, main=main_label)

lq_i = which(data[,"case_vs_control_de_lq"] == 1 & data[,"case_vs_control_de_hq"] == 0)
hq_i = which(data[,"case_vs_control_de_hq"] == 1)

points(x=log2(data[lq_i,"case_fpkm"]+variance_stabilization), y=log2(data[lq_i,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="orange")
points(x=log2(data[hq_i,"case_fpkm"]+variance_stabilization), y=log2(data[hq_i,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="red")
abline(a=0, b=1, lwd=2, lty=2, col="black")
abline(a=min_diff, b=1, lwd=1, lty=2, col="black")
abline(a=min_diff*-1, b=1, lwd=1, lty=2, col="black")
legend_text = c("no diff", "lq diff", "hq diff")
legend("topleft", legend_text, pch=16, bg="white", col=c("blue","orange","red"))
dev.off()

#Scatter plot of case vs. control FPKM values after normalization with hq and lq features indicated - protein_coding only
main_label = paste(case_label, " vs. ", control_label, " expression (protein coding ", type, ")", sep="")
x_label = paste(case_label, " - Cufflinks FPKM (log2)", sep="")
y_label = paste(control_label, " - Cufflinks FPKM (log2)", sep="")
png(filename = "case_vs_control_fpkm_scatter_postnorm_coding_de.png", width = 800, height = 800, bg = "white")
pc = which(data[,"biotype"] == "protein_coding")
plot(x=log2(data[pc,"case_fpkm"]+variance_stabilization), y=log2(data[pc,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="blue", xlab=x_label, ylab=y_label, main=main_label)

lq_i = which(data[,"case_vs_control_de_lq"] == 1 & data[,"case_vs_control_de_hq"] == 0 & data[,"biotype"] == "protein_coding")
hq_i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"biotype"] == "protein_coding")

points(x=log2(data[lq_i,"case_fpkm"]+variance_stabilization), y=log2(data[lq_i,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="orange")
points(x=log2(data[hq_i,"case_fpkm"]+variance_stabilization), y=log2(data[hq_i,"control_fpkm"]+variance_stabilization), pch=16, cex=0.75, col="red")
abline(a=0, b=1, lwd=2, lty=2, col="black")
abline(a=min_diff, b=1, lwd=1, lty=2, col="black")
abline(a=min_diff*-1, b=1, lwd=1, lty=2, col="black")
legend_text = c("no diff", "lq diff", "hq diff")
legend("topleft", legend_text, pch=16, bg="white", col=c("blue","orange","red"))
dev.off()


#Write out the output files

#- case_vs_control.tsv
write.table(data, file = "case_vs_control.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.lq.de.tsv
i = which(data[,"case_vs_control_de_lq"] == 1)
x = data[i,]
o = order(abs(x[,"case_vs_control_log2_de"]), decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.lq.de.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.lq.up.tsv
i = which(data[,"case_vs_control_de_lq"] == 1 & data[,"case_vs_control_log2_de"] > 0)
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.lq.up.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.lq.down.tsv
i = which(data[,"case_vs_control_de_lq"] == 1 & data[,"case_vs_control_log2_de"] < 0)
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=FALSE)
write.table(x[o,], file = "case_vs_control.lq.down.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.hq.de.tsv
i = which(data[,"case_vs_control_de_hq"] == 1)
x = data[i,]
o = order(abs(x[,"case_vs_control_log2_de"]), decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.hq.de.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.hq.up.tsv
i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"case_vs_control_log2_de"] > 0)
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.hq.up.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.hq.down.tsv
i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"case_vs_control_log2_de"] < 0)
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=FALSE)
write.table(x[o,], file = "case_vs_control.hq.down.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.coding.hq.de.tsv
i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"biotype"] == "protein_coding")
x = data[i,]
o = order(abs(x[,"case_vs_control_log2_de"]), decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.coding.hq.de.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.coding.hq.up.tsv
i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"case_vs_control_log2_de"] > 0 & data[,"biotype"] == "protein_coding")
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=TRUE)
write.table(x[o,], file = "case_vs_control.coding.hq.up.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#- case_vs_control.coding.hq.down.tsv
i = which(data[,"case_vs_control_de_hq"] == 1 & data[,"case_vs_control_log2_de"] < 0 & data[,"biotype"] == "protein_coding")
x = data[i,]
o = order(x[,"case_vs_control_log2_de"], decreasing=FALSE)
write.table(x[o,], file = "case_vs_control.coding.hq.down.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

