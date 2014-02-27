#!/usr/bin/env Rscript
#Written by Malachi Griffith
args = (commandArgs(TRUE))
working_dir = args[1] #Directory where output will be written
readcounts_file = args[2]  #snvs.hq.tier1.v1.annotated.compact.readcounts.tsv
gene_expression_file = args[3]  #isoforms.merged.fpkm.expsort.tsv

#Example execution
#WGS_vs_Exome_vs_RNAseq_VAF_and_FPKM.R  /gscmnt/sata132/techd/mgriffit/hgs/test/ /gscmnt/sata132/techd/mgriffit/hgs/all1/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv /gscmnt/sata132/techd/mgriffit/hgs/all1/rnaseq/tumor/absolute/isoforms_merged/isoforms.merged.fpkm.expsort.tsv

#working_dir = "/gscmnt/sata132/techd/mgriffit/luc/rnaseq_vs_snv/readcounts_merge_hg19/LUC1/summary/"
#readcounts_file = "/gscmnt/sata132/techd/mgriffit/hgs/all1/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv"
#gene_expression_file = "/gscmnt/sata132/techd/mgriffit/hgs/all1/rnaseq/tumor/absolute/isoforms_merged/isoforms.merged.fpkm.expsort.tsv"

#Load libraries
library(ggplot2)

if (length(args) < 2){
  message_text1 = "Required arguments missing for WGS_vs_Exome_vs_RNAseq_VAF_and_FPKM.R"
  stop(message_text1)
}

#Initialize the gene expression object, if the filepath was passed, load it, otherwise delete the object
gene_expression=NULL
if (is.na(args[3])){
  print ("Gene expression file not defined")
  rm(gene_expression)
}else{
  gene_expression=read.table(file=gene_expression_file, header=TRUE, sep="\t", as.is=c(1:4,10))
}

#Load the readcounts object
readcounts=read.table(file=readcounts_file, header=TRUE, sep="\t", as.is=c(1:6))

#Set the working directory
setwd(working_dir)

#Define some arbitrary cutoffs (TODO: be more systematic / rationalized)
var_stabilization = 0.1
min_support_cutoff = 1
rc_cutoff = 5
fpkm_cutoff = 1
vaf_diff_cutoff_pos = 20
vaf_diff_cutoff_neg = vaf_diff_cutoff_pos*-1
psize = 3

#Create a series of X-Y scatter plots of Tumor/Normal Variant Allele Frequencies (VAF) using ggplot2
#Plot WGS vs. RNA-seq, Exome vs. RNA-seq, and WGS vs. Exome according to data availability
#If available, plot RNA-seq Tumor VAF vs. RNA-seq Normal VAF

#For plots involving RNAseq VAFs ... color points on a yellow->red spectrum according to Gene FPKM
#If NA's are present in these data (genes where the mutation could not be mapped successfully to their gene FPKMs from the RNA-seq analysis) ...
#Reset these data points to 0.  This is not ideal!  Would be better to have complete mapping of genes!
#If some non-NA value are present, only these will be used for the plots
if (length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  i = which(is.na(readcounts[,"RNAseq_Tumor_gene_FPKM"]))
  if (length(i)){
    readcounts[i,"RNAseq_Tumor_gene_FPKM"] = 0
  }
}
if (length(which(names(readcounts)=="RNAseq_Normal_gene_FPKM"))){
  i = which(is.na(readcounts[,"RNAseq_Normal_gene_FPKM"]))
  if (length(i)){
    readcounts[i,"RNAseq_Normal_gene_FPKM"] = 0
  }
}


########################################################################################################
ypos=115
#1.) Tumor VAF - WGS vs RNA-seq - Check for existence of these data columns
if (length(which(names(readcounts)=="WGS_Tumor_VAF")) & length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  #A.) All data
  #Adjust the Tumor Gene FPKM values
  Gene_FPKM=log2(readcounts[,"RNAseq_Tumor_gene_FPKM"] + var_stabilization)

  if (min(Gene_FPKM, na.rm=TRUE) < 0){
	  Gene_FPKM = Gene_FPKM + abs(min(Gene_FPKM, na.rm=TRUE))
  }
  n = length(readcounts[,"WGS_Tumor_VAF"])
  m = length(which(readcounts[,"WGS_Tumor_VAF"] > 0 & readcounts[,"RNAseq_Tumor_VAF"] > 0))
  p = round(((m/n)*100), digits=1)
  spearmanR = cor(x=readcounts[,"WGS_Tumor_VAF"], y=readcounts[,"RNAseq_Tumor_VAF"], method="spearman")
  pearsonR = cor(x=readcounts[,"WGS_Tumor_VAF"], y=readcounts[,"RNAseq_Tumor_VAF"], method="pearson")
  text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
  text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
  text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")
  
  pdf(file="Tumor_VAF_WGS_vs_RNAseq_scatter.pdf")
  d = ggplot(data=readcounts, aes(x=WGS_Tumor_VAF, y=RNAseq_Tumor_VAF))
  d = d + geom_point(aes(colour=Gene_FPKM), size=psize) + scale_colour_gradient('FPKM', low="yellow", high="red")
  d = d + xlab("Tumor Variant Allele Frequency (WGS)") + ylab("Tumor Variant Allele Frequency (RNA-seq)")
  d = d + xlim(c(0,100)) + ylim(c(0,ypos))
  if (packageVersion("ggplot2") <= "0.8.9"){
    d = d + opts(title="Detection and RNA expression of variants identified by WGS")
  }else{
    d = d + labs(title="Detection and RNA expression of variants identified by WGS")
  }
  d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
  d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
  d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1)
  d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))
  print(d)
  dev.off()
  
  #B.) Read count filtered data
  #Apply a minimum total tumor read count cutoff (ref + var) and then replot the data - WGS vs RNA-seq
  #- RNAseq and WGS Read counts covering the variant position must be higher than a min read cutoff
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) > rc_cutoff) & ((readcounts[,"WGS_Tumor_ref_rc"] + readcounts[,"WGS_Tumor_var_rc"]) > rc_cutoff))
  if (length(i)){
    readcounts_cut = readcounts[i,]
    Gene_FPKM=log2(readcounts_cut[,"RNAseq_Tumor_gene_FPKM"] + var_stabilization)
    if (min(Gene_FPKM, na.rm=TRUE) < 0){
	    Gene_FPKM = Gene_FPKM + abs(min(Gene_FPKM, na.rm=TRUE))
    }
    n = length(readcounts_cut[,"WGS_Tumor_VAF"])
    m = length(which(readcounts_cut[,"WGS_Tumor_VAF"] > 0 & readcounts_cut[,"RNAseq_Tumor_VAF"] > 0))
    p = round(((m/n)*100), digits=1)
    spearmanR = cor(x=readcounts_cut[,"WGS_Tumor_VAF"], y=readcounts_cut[,"RNAseq_Tumor_VAF"], method="spearman")
    pearsonR = cor(x=readcounts_cut[,"WGS_Tumor_VAF"], y=readcounts_cut[,"RNAseq_Tumor_VAF"], method="pearson")
    text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
    text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
    text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")

    pdf(file="Tumor_VAF_WGS_vs_RNAseq_ReadCutoff_scatter.pdf")
    d = ggplot(data=readcounts_cut, aes(x=WGS_Tumor_VAF, y=RNAseq_Tumor_VAF)) 
    d = d + geom_point(aes(colour=Gene_FPKM), size=psize) + scale_colour_gradient('FPKM', low="yellow", high="red") 
    d = d + xlab("Tumor Variant Allele Frequency (WGS)") + ylab("Tumor Variant Allele Frequency (RNA-seq)")
    d = d + xlim(c(0,100)) + ylim(c(0,ypos))
    if (packageVersion("ggplot2") <= "0.8.9"){
      d = d + opts(title="Detection and RNA expression of variants identified by WGS (> rc)")
    }else{
      d = d + labs(title="Detection and RNA expression of variants identified by WGS (> rc)")      
    }
    d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
    d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
    d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1)
	  d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))
    print(d)
    dev.off()
  }
}


########################################################################################################
#2.) Tumor VAF - Exome vs RNA-seq - Check for existence of these data columns
if (length(which(names(readcounts)=="Exome_Tumor_VAF")) & length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  #A.) All data
  #Adjust the Tumor Gene FPKM values
  Gene_FPKM=log2(readcounts[,"RNAseq_Tumor_gene_FPKM"] + var_stabilization)
    if (min(Gene_FPKM, na.rm=TRUE) < 0){
	  Gene_FPKM = Gene_FPKM + abs(min(Gene_FPKM, na.rm=TRUE))
  }
  n = length(readcounts[,"Exome_Tumor_VAF"])
  m = length(which(readcounts[,"Exome_Tumor_VAF"] > 0 & readcounts[,"RNAseq_Tumor_VAF"] > 0))
  p = round(((m/n)*100), digits=1)
  spearmanR = cor(x=readcounts[,"Exome_Tumor_VAF"], y=readcounts[,"RNAseq_Tumor_VAF"], method="spearman")
  pearsonR = cor(x=readcounts[,"Exome_Tumor_VAF"], y=readcounts[,"RNAseq_Tumor_VAF"], method="pearson")
  text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
  text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
  text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")

  pdf(file="Tumor_VAF_Exome_vs_RNAseq_scatter.pdf")
  d = ggplot(data=readcounts, aes(x=Exome_Tumor_VAF, y=RNAseq_Tumor_VAF))
  d = d + geom_point(aes(colour=Gene_FPKM), size=psize) + scale_colour_gradient('FPKM', low="yellow", high="red")
  d = d + xlab("Tumor Variant Allele Frequency (Exome)") + ylab("Tumor Variant Allele Frequency (RNA-seq)")
  d = d + xlim(c(0,100)) + ylim(c(0,ypos))
  if (packageVersion("ggplot2") <= "0.8.9"){
    d = d + opts(title="Detection and RNA expression of variants identified by Exome")
  }else{
    d = d + labs(title="Detection and RNA expression of variants identified by Exome")
  }
  d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
  d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
  d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1)
	d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))
  print(d)
  dev.off()
  
  #B.) Read count filtered data
  #Apply a minimum total tumor read count cutoff (ref + var) and then replot the data - Exome vs RNA-seq
  #- RNAseq and Exome Read counts covering the variant position must be higher than a min read cutoff
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) > rc_cutoff) & ((readcounts[,"Exome_Tumor_ref_rc"] + readcounts[,"Exome_Tumor_var_rc"]) > rc_cutoff))
  if (length(i)){
    readcounts_cut = readcounts[i,]
    Gene_FPKM=log2(readcounts_cut[,"RNAseq_Tumor_gene_FPKM"] + var_stabilization)
    if (min(Gene_FPKM, na.rm=TRUE) < 0){
	    Gene_FPKM = Gene_FPKM + abs(min(Gene_FPKM, na.rm=TRUE))
    }
    n = length(readcounts_cut[,"Exome_Tumor_VAF"])
    m = length(which(readcounts_cut[,"Exome_Tumor_VAF"] > 0 & readcounts_cut[,"RNAseq_Tumor_VAF"] > 0))
    p = round(((m/n)*100), digits=1)
    spearmanR = cor(x=readcounts_cut[,"Exome_Tumor_VAF"], y=readcounts_cut[,"RNAseq_Tumor_VAF"], method="spearman")
    pearsonR = cor(x=readcounts_cut[,"Exome_Tumor_VAF"], y=readcounts_cut[,"RNAseq_Tumor_VAF"], method="pearson")
    text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
    text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
    text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")

    pdf(file="Tumor_VAF_Exome_vs_RNAseq_ReadCutoff_scatter.pdf")
    d = ggplot(data=readcounts_cut, aes(x=Exome_Tumor_VAF, y=RNAseq_Tumor_VAF)) 
    d = d + geom_point(aes(colour=Gene_FPKM), size=psize) + scale_colour_gradient('FPKM', low="yellow", high="red")
    d = d + xlab("Tumor Variant Allele Frequency (Exome)") + ylab("Tumor Variant Allele Frequency (RNA-seq)")
    d = d + xlim(c(0,100)) + ylim(c(0,ypos))
    if (packageVersion("ggplot2") <= "0.8.9"){
      d = d + opts(title="Detection and RNA expression of variants identified by Exome (> rc)")
    }else{
      d = d + labs(title="Detection and RNA expression of variants identified by Exome (> rc)")
    }
    d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
    d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
    d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1)
	  d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))    
    print(d)
    dev.off()
  }
}


########################################################################################################
#3.) Tumor VAF - WGS vs Exome - Check for existence of these data columns
if (length(which(names(readcounts)=="WGS_Tumor_VAF")) & length(which(names(readcounts)=="Exome_Tumor_VAF"))){
  #A.) All data

  n = length(readcounts[,"WGS_Tumor_VAF"])
  m = length(which(readcounts[,"WGS_Tumor_VAF"] > 0 & readcounts[,"Exome_Tumor_VAF"] > 0))
  p = round(((m/n)*100), digits=1)
  spearmanR = cor(x=readcounts[,"WGS_Tumor_VAF"], y=readcounts[,"Exome_Tumor_VAF"], method="spearman")
  pearsonR = cor(x=readcounts[,"WGS_Tumor_VAF"], y=readcounts[,"Exome_Tumor_VAF"], method="pearson")
  text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
  text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
  text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")

  pdf(file="Tumor_VAF_WGS_vs_Exome_scatter.pdf")
  d = ggplot(data=readcounts, aes(x=WGS_Tumor_VAF, y=Exome_Tumor_VAF))
  d = d + geom_point(size=psize)
  d = d + xlab("Tumor Variant Allele Frequency (WGS)") + ylab("Tumor Variant Allele Frequency (Exome)")
  d = d + xlim(c(0,100)) + ylim(c(0,ypos))
  if (packageVersion("ggplot2") <= "0.8.9"){
    d = d + opts(title="Detection of variants by WGS vs. Exome")
  }else{
    d = d + labs(title="Detection of variants by WGS vs. Exome")
  }
  d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
  d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
  d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1) 
	d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))
  print(d)
  dev.off()
  
  #B.) Read count filtered data
  #Apply a minimum total tumor read count cutoff (ref + var) and then replot the data - WGS vs Exome
  #- Exome and WGS Read counts covering the variant position must be higher than a min read cutoff
  i = which(((readcounts[,"Exome_Tumor_ref_rc"] + readcounts[,"Exome_Tumor_var_rc"]) > rc_cutoff) & ((readcounts[,"WGS_Tumor_ref_rc"] + readcounts[,"WGS_Tumor_var_rc"]) > rc_cutoff))
  if (length(i)){
    readcounts_cut = readcounts[i,]

    n = length(readcounts_cut[,"WGS_Tumor_VAF"])
    m = length(which(readcounts_cut[,"WGS_Tumor_VAF"] > 0 & readcounts_cut[,"Exome_Tumor_VAF"] > 0))
    p = round(((m/n)*100), digits=1)
    spearmanR = cor(x=readcounts_cut[,"WGS_Tumor_VAF"], y=readcounts_cut[,"Exome_Tumor_VAF"], method="spearman")
    pearsonR = cor(x=readcounts_cut[,"WGS_Tumor_VAF"], y=readcounts_cut[,"Exome_Tumor_VAF"], method="pearson")
    text1 = paste("n = ", n, " (", m, " observed in both", " [", p, "%])", sep="")
    text2 = paste("R = ", round(spearmanR, digits=2), " (spearman)", sep="")
    text3 = paste("R = ", round(pearsonR, digits=2), " (pearson)", sep="")

    pdf(file="Tumor_VAF_WGS_vs_Exome_ReadCutoff_scatter.pdf")
    d = ggplot(data=readcounts_cut, aes(x=WGS_Tumor_VAF, y=Exome_Tumor_VAF)) 
    d = d + geom_point(size=psize)
    d = d + xlab("Tumor Variant Allele Frequency (WGS)") + ylab("Tumor Variant Allele Frequency (Exome)")
    d = d + xlim(c(0,100)) + ylim(c(0,ypos))
    if (packageVersion("ggplot2") <= "0.8.9"){
      d = d + opts(title="Detection of variants by WGS vs. Exome (> rc)")
    }else{
      d = d + labs(title="Detection of variants by WGS vs. Exome (> rc)")
    }
    d = d + geom_abline(intercept=0, slope=1, linetype=2, size=0.2)
    d = d + geom_abline(intercept=20, slope=1, linetype=2, size=0.1)
    d = d + geom_abline(intercept=-20, slope=1, linetype=2, size=0.1) 
	  d = d + geom_text(aes(x2,y2,label = textlist, hjust=0), data.frame(x2=c(0,0,0), y2=c(ypos,ypos-5,ypos-10), textlist=c(text1,text2,text3)))
    print(d)
    dev.off()
  }
}

#Plot distributions of Tumor FPKM values for Mutant vs. Non-Mutant genes
if (exists("gene_expression")){
  if (length(which(names(gene_expression)=="FPKM"))){
  
    #Transform the FPKM values for display purposes
	  gene_expression[,"FPKM_trans"] = log2(gene_expression[,"FPKM"] + var_stabilization)

    #Identify the subset of all genes that are mutant
    mutant_gene_names = unique(readcounts[,"mapped_gene_name"])
    gene_expression[,"Group"] = "Non-Mutant"
    i = which(gene_expression[,"mapped_gene_name"] %in% mutant_gene_names)
    gene_expression[i,"Group"] = "Mutant"
    mutant_count = length(i)
    title_text = paste("Expression of non-mutant vs. mutant genes (n = ", mutant_count, ")", sep="")

    #Plot distribution of all Gene FPKMs vs. those of mutant genes - Tumor
    #Boxplot
    pdf("Tumor_MutantGene_vs_NonMutantGene_FPKM_boxplot.pdf")
    p = ggplot(gene_expression, aes(factor(Group), FPKM_trans)) 
    p = p + geom_boxplot(aes(fill=Group)) + xlab("Gene group") + ylab("FPKM (log2)")
    if (packageVersion("ggplot2") <= "0.8.9"){
      p = p + opts(title=title_text)
    }else{
      p = p + labs(title=title_text)
    }
    print(p)
    dev.off()

    #Density plot
    pdf("Tumor_MutantGene_vs_NonMutantGene_FPKM_density.pdf")
    p = ggplot(gene_expression, aes(FPKM_trans, fill = Group)) + geom_density(alpha = 0.2) + xlab("FPKM (log2)") + ylab("Density")
    if (packageVersion("ggplot2") <= "0.8.9"){
      p = p + opts(title=title_text)
    }else{
      p = p + labs(title=title_text)
    }
    print(p)
    dev.off()

    #Repeat these plots using only genes with nonzero expression
    exp_i = which(gene_expression[,"FPKM"] > 0)
    ge_i = gene_expression[exp_i,]
    mutant_count = length(which(ge_i[,"Group"] == "Mutant"))
    title_text = paste("Expression of expressed non-mutant vs. mutant genes (n = ", mutant_count, ")", sep="")

    pdf("Tumor_MutantGene_vs_NonMutantGene_FPKM_ReadCutoff_boxplot.pdf")
    p = ggplot(ge_i, aes(factor(Group), FPKM_trans)) 
    p = p + geom_boxplot(aes(fill=Group)) + xlab("Gene group") + ylab("FPKM (log2) (genes with FPKM > 0)")
    if (packageVersion("ggplot2") <= "0.8.9"){
      p = p + opts(title=title_text)
    }else{
      p = p + labs(title=title_text)
    }
    print(p)
    dev.off()

    pdf("Tumor_MutantGene_vs_NonMutantGene_FPKM_ReadCutoff_density.pdf")
    p = ggplot(ge_i, aes(FPKM_trans, fill = Group)) + geom_density(alpha = 0.2) + xlab("FPKM (log2) (genes with FPKM > 0)") + ylab("Density") 
    if (packageVersion("ggplot2") <= "0.8.9"){
      p = p + opts(title=title_text)
    }else{
      p = p + labs(title=title_text)
    }
    print(p)
    dev.off()
  }
}

#Plot VAF values as histograms
plotVafHist = function(dataset, sample_type, data_type){
  colname = paste(data_type, "_", sample_type, "_VAF", sep="")
  filename = paste(sample_type, "_VAF_", data_type, "_hist.pdf", sep="")
  title_text = paste("Distribution of variant allele frequencies (", sample_type, " tissue, ", data_type, ")", sep="")
  x_label = paste(data_type, " ", sample_type, " variant allele frequency", sep="")
  if (length(which(names(dataset)==colname))){
    y_label = paste("Density (n = ", length(dataset[,colname]), " variants)", sep="")
    dataset[,"VAF"] = dataset[,colname]
    pdf(filename)
    m = ggplot(dataset, aes(x=VAF))
    m = m + geom_histogram(aes(y = ..density.., fill= ..count..))
    m = m + geom_density()
    if (packageVersion("ggplot2") <= "0.8.9"){
      m = m + opts(title=title_text) 
    }else{
      m = m + labs(title=title_text) 
    }
    m = m + xlab(x_label) + ylab(y_label)
    print(m)
    dev.off()
  }
}
plotVafHist(readcounts, "Normal", "WGS")
plotVafHist(readcounts, "Tumor", "WGS")
plotVafHist(readcounts, "Normal", "Exome")
plotVafHist(readcounts, "Tumor", "Exome")
plotVafHist(readcounts, "Normal", "RNAseq")
plotVafHist(readcounts, "Tumor", "RNAseq")


#Create density plots for the distributions of Tumor VAFs (and Normal?) for:
#WGS + Exome + RNAseq (all on one plot)
x=NULL
y=NULL
z=NULL
var_count = 0
if (length(which(names(readcounts)=="WGS_Tumor_VAF"))){
  x=readcounts[,"WGS_Tumor_VAF"]
  var_count = length(x)
}
if (length(which(names(readcounts)=="Exome_Tumor_VAF"))){
  y=readcounts[,"Exome_Tumor_VAF"]
  var_count = length(y)
}
if (length(which(names(readcounts)=="RNAseq_Tumor_VAF"))){
  z=readcounts[,"RNAseq_Tumor_VAF"]
  var_count = length(z)
}
classes = c(rep("WGS", length(x)), rep("Exome", length(y)), rep("RNAseq", length(z)))
vafs = data.frame(c(x,y,z), classes)
names(vafs) = c("Tumor_VAF","Group")

#Only plot if at least one category has at least 3 values
min_cat = min(table(vafs[,"Group"]))
if (min_cat >= 3){
  y_label = paste("Density (n = ", var_count, " variants)", sep="")
  pdf("Tumor_VAF_AllDataSources_density.pdf")
  p = ggplot(vafs, aes(Tumor_VAF, fill = Group)) + geom_density(alpha = 0.2) + xlab("Tumor Variant Allele Frequency") + ylab(y_label)
  if (packageVersion("ggplot2") <= "0.8.9"){
    p = p + opts(title="Comparison of variant allele frequencies from each data source")
  }else{
    p = p + labs(title="Comparison of variant allele frequencies from each data source")
  }
  print(p)
  dev.off()
}

#Summarize coverage levels for mutation positions according to: WGS, Exome, and RNA-seq 
plotCoverageHist = function(data, ref_rc, var_rc, filename, title_name){
  #With outliers
  filename1 = paste(filename, "_hist.pdf", sep="")
  z = data[,c(ref_rc,var_rc)]
  z[,"coverage"] = z[,ref_rc]+z[,var_rc]
  median_cov = median(z[,"coverage"], na.rm=TRUE)
  x_label = paste("Read coverage (median X = ", round(median_cov, digits=1), ")", sep="")
  y_label = paste("Density (n = ", dim(z)[1], " variant positions)", sep="")
  pdf(filename1)
  m = ggplot(z, aes(x=coverage))
  m = m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density()
  if (packageVersion("ggplot2") <= "0.8.9"){
    m = m + opts(title=title_name) 
  }else{
    m = m + labs(title=title_name) 
  }
  m = m + xlab(x_label) + ylab(y_label)
  print(m)
  dev.off()

  #Remove outliers
  outliers=boxplot(z[,"coverage"], plot=FALSE)$out
  z = z[which(!z[,"coverage"] %in% outliers),]
  filename2 = paste(filename, "_RmOutliers_hist.pdf", sep="")
  x_label = paste("Read coverage (median X = ", round(median_cov, digits=1), ") - Outliers removed", sep="")
  y_label = paste("Density (n = ", dim(z)[1], " variant positions)", sep="")
  pdf(filename2)
  m = ggplot(z, aes(x=coverage))
  m = m + geom_histogram(aes(y = ..density.., fill= ..count..)) + geom_density()
  if (packageVersion("ggplot2") <= "0.8.9"){
    m = m + opts(title=title_name) 
  }else{
    m = m + labs(title=title_name) 
  }
  m = m + xlab(x_label) + ylab(y_label)
  print(m)
  dev.off()
}
if (length(which(names(readcounts)=="WGS_Tumor_ref_rc"))){
  plotCoverageHist(readcounts, "WGS_Tumor_ref_rc", "WGS_Tumor_var_rc", "Tumor_VariantReadCoverage_WGS", "Distribution of variant read coverage (Tumor, WGS)")
}
if (length(which(names(readcounts)=="WGS_Normal_ref_rc"))){
  plotCoverageHist(readcounts, "WGS_Normal_ref_rc", "WGS_Normal_var_rc", "Normal_VariantReadCoverage_WGS", "Distribution of variant read coverage (Normal, WGS)")
}
if (length(which(names(readcounts)=="Exome_Tumor_ref_rc"))){
  plotCoverageHist(readcounts, "Exome_Tumor_ref_rc", "Exome_Tumor_var_rc", "Tumor_VariantReadCoverage_Exome", "Distribution of variant read coverage (Tumor, Exome)")
}
if (length(which(names(readcounts)=="Exome_Normal_ref_rc"))){
  plotCoverageHist(readcounts, "Exome_Normal_ref_rc", "Exome_Normal_var_rc", "Normal_VariantReadCoverage_Exome", "Distribution of variant read coverage (Normal, Exome)")
}
if (length(which(names(readcounts)=="RNAseq_Tumor_ref_rc"))){
  plotCoverageHist(readcounts, "RNAseq_Tumor_ref_rc", "RNAseq_Tumor_var_rc", "Tumor_VariantReadCoverage_RNAseq", "Distribution of variant read coverage (Tumor, RNAseq)")
}
if (length(which(names(readcounts)=="RNAseq_Normal_ref_rc"))){
  plotCoverageHist(readcounts, "RNAseq_Normal_ref_rc", "RNAseq_Normal_var_rc", "Normal_VariantReadCoverage_RNAseq", "Distribution of variant read coverage (Normal, RNAseq)")
}

#Same idea as above but combine on one plot
#Tumor
plotCoverageComparison = function(dataset, sample_type, remove_outliers){
  x=NULL; xt=NULL; y=NULL; yt=NULL; z=NULL; zt=NULL; xmedian=NULL; ymedian=NULL; zmedian=NULL;
  wgs_ref_rc = paste("WGS_", sample_type, "_ref_rc", sep="")
  wgs_var_rc = paste("WGS_", sample_type, "_var_rc", sep="")
  exome_ref_rc = paste("Exome_", sample_type, "_ref_rc", sep="")
  exome_var_rc = paste("Exome_", sample_type, "_var_rc", sep="")
  rnaseq_ref_rc = paste("RNAseq_", sample_type, "_ref_rc", sep="")
  rnaseq_var_rc = paste("RNAseq_", sample_type, "_var_rc", sep="")
  var_count = 0
  if (length(which(names(dataset)==wgs_ref_rc))){
    x=dataset[,wgs_ref_rc]+dataset[,wgs_var_rc]
    xmedian=median(x)
    if (remove_outliers){
      outliers=boxplot(x, plot=FALSE)$out
      x = x[which(!x %in% outliers)]
    }
    xt=paste("WGS=", round(xmedian, digits=1), "X ", sep="")
    if (length(x) > var_count){var_count = length(x)}
  }
  if (length(which(names(dataset)==exome_ref_rc))){
    y=dataset[,exome_ref_rc]+dataset[,exome_var_rc]
    ymedian=median(y)
    if (remove_outliers){
      outliers=boxplot(y, plot=FALSE)$out
      y = y[which(!y %in% outliers)]
    }
    yt=paste("Exome=" , round(ymedian, digits=1), "X ", sep="")
    if (length(y) > var_count){var_count = length(y)}
  }
  #For RNAseq only, calculate the median from nonzero values only
  if (length(which(names(dataset)==rnaseq_ref_rc))){
    z=dataset[,rnaseq_ref_rc]+dataset[,rnaseq_var_rc]
    z_nonzero = z[which(z > 0)]
    zmedian=median(z_nonzero)
    if (remove_outliers){
      outliers=boxplot(z, plot=FALSE)$out
      z = z[which(!z %in% outliers)]
    }
    zt=paste("RNAseq=", round(zmedian, digits=1), "X ", sep="")
    if (length(z) > var_count){var_count = length(z)}
  }
  classes = c(rep("WGS", length(x)), rep("Exome", length(y)), rep("RNAseq", length(z)))
  rcs = data.frame(c(x,y,z), classes)
  names(rcs) = c("ReadCounts","DataType")

  min_cat = min(table(rcs[,"DataType"]))
  if (min_cat >= 3){
    x_label = paste(sample_type, " Read Coverage ( ", zt, xt, yt, ")", sep="")
    y_label = paste("Density (n = ", var_count, " variants)", sep="")
    filename = paste(sample_type,"_ReadCoverage_AllDataSources_density.pdf", sep="")
    if (remove_outliers){
   	  filename = paste(sample_type,"_ReadCoverage_AllDataSources_RmOutliers_density.pdf", sep="")	
    }
    pdf(filename)
    
    p = ggplot(rcs, aes(ReadCounts, fill = DataType)) + geom_density(alpha = 0.2)
    p = p + xlab(x_label) + ylab(y_label)
    if (packageVersion("ggplot2") <= "0.8.9"){
      p = p + opts(title="Comparison of read coverages from each data source")
    }else{
      p = p + labs(title="Comparison of read coverages from each data source")
    }
    p = p + geom_vline(xintercept=xmedian, linetype=3, color="black", size=0.3) + geom_vline(xintercept=ymedian, linetype=3, color="black", size=0.3) + geom_vline(xintercept=zmedian, linetype=3, color="black", size=0.3)
    print(p)
    dev.off()
  }
}
plotCoverageComparison(readcounts, "Tumor", remove_outliers=FALSE)
plotCoverageComparison(readcounts, "Normal", remove_outliers=FALSE)
plotCoverageComparison(readcounts, "Tumor", remove_outliers=TRUE)
plotCoverageComparison(readcounts, "Normal", remove_outliers=TRUE)

#Plot expression values from the tumor - FPKMs vs. Coverage values (and then mark mutant genes)
if (exists("gene_expression")){
  if (length(which(names(gene_expression)=="FPKM"))){
    #Transform the FPKM values for display purposes
	gene_expression[,"FPKM_trans"] = log2(gene_expression[,"FPKM"] + var_stabilization)
	#Similar transformation on the coverage values
	gene_expression[,"coverage_trans"] = log2(gene_expression[,"coverage"] + var_stabilization)
	#Identify the subset of all genes that are mutant
	mutant_gene_names = unique(readcounts[,"mapped_gene_name"])
	gene_expression[,"Group"] = "Non-Mutant"
	i = which(gene_expression[,"mapped_gene_name"] %in% mutant_gene_names)
	gene_expression[i,"Group"] = "Mutant"
	mutant_count = length(i)
	pointsize=1
	if (mutant_count > 150){
      pointsize=0.35	
	}
    jpeg("Tumor_Gene_Expression_Mutant_vs_NonMutant_scatter.jpg", quality=100, width=800, height=800)
	plot(x=gene_expression[,"FPKM_trans"], y=gene_expression[,"coverage_trans"], col="blue", pch=16, cex=0.2, 
	     xlab="FPKM (log2)", ylab="Coverage (log2)", main="Gene expression value for mutant and non-mutant genes")
    points(x=gene_expression[i,"FPKM_trans"], y=gene_expression[i,"coverage_trans"], col="magenta", pch=16, cex=1)
	legend("topleft", c("Mutant", "Non-Mutant"), col=c("magenta", "blue"), pch=16)
    dev.off()
  }
}


#SUMMARY STATISTICS
#Write out all stats in the following format:
#Question | Answer | Data type | analysis type | statistic type | Extra description

#Initialize a stats data.frame
stats = data.frame(NA,NA,NA,NA,NA,NA)
names(stats) = c("Question", "Answer", "Data_Type", "Analysis_Type", "Statistic_Type", "Extra_Description")
total_variants = dim(readcounts)[1]
stats[dim(stats)[1],] = c("Number of variants", total_variants,"RNA-seq","SNV","Count","Variants from WGS/Exome")

#Correlation values for:
# - WGS vs. Exome VAFs
if (length(which(names(readcounts)=="Exome_Tumor_VAF")) & length(which(names(readcounts)=="WGS_Tumor_VAF"))){
  i = which(readcounts[,"Exome_Tumor_VAF"] > rc_cutoff & readcounts[,"WGS_Tumor_VAF"] > rc_cutoff)
  wgs_vs_exome_vaf_cor = cor(x=readcounts[,"Exome_Tumor_VAF"], y=readcounts[,"WGS_Tumor_VAF"], method="spearman")
  wgs_vs_exome_vaf_cor_rc = cor(x=readcounts[i,"Exome_Tumor_VAF"], y=readcounts[i,"WGS_Tumor_VAF"], method="spearman")
  stats[dim(stats)[1]+1,] = c("Exome vs. WGS Tumor VAF correlation", wgs_vs_exome_vaf_cor, "WGS and Exome","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for Exome vs. WGS data")
  stats[dim(stats)[1]+1,] = c("Exome vs. WGS Tumor VAF correlation - with min coverage cutoff", wgs_vs_exome_vaf_cor_rc, "WGS and Exome","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for Exome vs. WGS data")
}
# - WGS vs. RNAseq VAFs (all & 'expressed' variants only)
if (length(which(names(readcounts)=="RNAseq_Tumor_VAF")) & length(which(names(readcounts)=="WGS_Tumor_VAF"))){
  i = which(readcounts[,"RNAseq_Tumor_VAF"] > rc_cutoff & readcounts[,"WGS_Tumor_VAF"] > rc_cutoff)
  rnaseq_vs_wgs_vaf_cor = cor(x=readcounts[,"RNAseq_Tumor_VAF"], y=readcounts[,"WGS_Tumor_VAF"], method="spearman")
  rnaseq_vs_wgs_vaf_cor_rc = cor(x=readcounts[i,"RNAseq_Tumor_VAF"], y=readcounts[i,"WGS_Tumor_VAF"], method="spearman")
  stats[dim(stats)[1]+1,] = c("RNAseq vs. WGS Tumor VAF correlation", rnaseq_vs_wgs_vaf_cor, "RNAseq and WGS","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for RNAseq vs. WGS data")
  stats[dim(stats)[1]+1,] = c("RNAseq vs. WGS Tumor VAF correlation - with min coverage cutoff", rnaseq_vs_wgs_vaf_cor_rc, "RNAseq and WGS","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for RNAseq vs. WGS data")
}
# - Exome vs. RNAseq VAFs (all & 'expressed' variants only)
if (length(which(names(readcounts)=="RNAseq_Tumor_VAF")) & length(which(names(readcounts)=="Exome_Tumor_VAF"))){
  i = which(readcounts[,"RNAseq_Tumor_VAF"] > rc_cutoff & readcounts[,"Exome_Tumor_VAF"] > rc_cutoff)
  exome_vs_rnaseq_vaf_cor = cor(x=readcounts[,"Exome_Tumor_VAF"], y=readcounts[,"RNAseq_Tumor_VAF"], method="spearman")
  exome_vs_rnaseq_vaf_cor_rc = cor(x=readcounts[i,"Exome_Tumor_VAF"], y=readcounts[i,"RNAseq_Tumor_VAF"], method="spearman")
  stats[dim(stats)[1]+1,] = c("RNAseq vs. Exome Tumor VAF correlation", exome_vs_rnaseq_vaf_cor, "RNAseq and Exome","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for RNAseq vs. Exome data")
  stats[dim(stats)[1]+1,] = c("RNAseq vs. Exome Tumor VAF correlation - with min coverage cutoff", exome_vs_rnaseq_vaf_cor_rc, "RNAseq and Exome","SNV","Spearman correlation", "Correlation between tumor variant allele frequencies for RNAseq vs. Exome data") 
}

# - Median coverage of variant positions in 
# - WGS Tumor
if (length(which(names(readcounts)=="WGS_Tumor_var_rc"))){
  wgs_median_tumor_coverage = median(readcounts[,"WGS_Tumor_ref_rc"] + readcounts[,"WGS_Tumor_var_rc"])
  wgs_median_tumor_vaf = median(readcounts[,"WGS_Tumor_VAF"])
  stats[dim(stats)[1]+1,] = c("WGS median tumor read coverage", wgs_median_tumor_coverage, "WGS","SNV", "Median", "Median tumor coverage (reads supporting reference and variant bases)") 
  stats[dim(stats)[1]+1,] = c("WGS median tumor VAF", wgs_median_tumor_vaf, "WGS","SNV", "Median", "Median tumor variant allele frequency") 
}
# - WGS Normal
if (length(which(names(readcounts)=="WGS_Normal_var_rc"))){
  wgs_median_normal_coverage = median(readcounts[,"WGS_Normal_ref_rc"] + readcounts[,"WGS_Normal_var_rc"])
  wgs_median_normal_vaf = median(readcounts[,"WGS_Normal_VAF"])
  stats[dim(stats)[1]+1,] = c("WGS median normal read coverage", wgs_median_normal_coverage, "WGS","SNV", "Median", "Median normal coverage (reads supporting reference and variant bases)")
  stats[dim(stats)[1]+1,] = c("WGS median normal VAF", wgs_median_normal_vaf, "WGS","SNV", "Median", "Median normal variant allele frequency")
}
# - Exome Tumor
if (length(which(names(readcounts)=="Exome_Tumor_var_rc"))){
  exome_median_tumor_coverage = median(readcounts[,"Exome_Tumor_ref_rc"] + readcounts[,"Exome_Tumor_var_rc"])
  exome_median_tumor_vaf = median(readcounts[,"Exome_Tumor_VAF"])
  stats[dim(stats)[1]+1,] = c("Exome median tumor read coverage", exome_median_tumor_coverage, "Exome","SNV", "Median", "Median tumor coverage (reads supporting reference and variant bases)") 
  stats[dim(stats)[1]+1,] = c("Exome median tumor VAF", exome_median_tumor_vaf, "Exome","SNV", "Median", "Median tumor variant allele frequency") 
}
# - Exome Normal
if (length(which(names(readcounts)=="Exome_Normal_var_rc"))){
  exome_median_normal_coverage = median(readcounts[,"Exome_Normal_ref_rc"] + readcounts[,"Exome_Normal_var_rc"])
  exome_median_normal_vaf = median(readcounts[,"Exome_Normal_VAF"])
  stats[dim(stats)[1]+1,] = c("Exome median normal read coverage", exome_median_normal_coverage, "Exome","SNV", "Median", "Median normal coverage (reads supporting reference and variant bases)") 
  stats[dim(stats)[1]+1,] = c("Exome median normal VAF", exome_median_normal_vaf, "Exome","SNV", "Median", "Median normal variant allele frequency") 
}

# - RNAseq Tumor (all, variant observed only, expressed gene only)
if (length(which(names(readcounts)=="RNAseq_Tumor_var_rc")) & length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  rnaseq_median_tumor_coverage = median(readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"])
  rnaseq_median_tumor_vaf = median(readcounts[,"RNAseq_Tumor_VAF"])
  i = which(readcounts[,"RNAseq_Tumor_var_rc"] > 0)
  rnaseq_median_tumor_coverage_obs = median(readcounts[i,"RNAseq_Tumor_ref_rc"] + readcounts[i,"RNAseq_Tumor_var_rc"])
  rnaseq_median_tumor_vaf_obs = median(readcounts[i,"RNAseq_Tumor_VAF"])
  i = which(readcounts[,"RNAseq_Tumor_gene_FPKM"] > fpkm_cutoff)
  rnaseq_median_tumor_coverage_exp = median(readcounts[i,"RNAseq_Tumor_ref_rc"] + readcounts[i,"RNAseq_Tumor_var_rc"])
  rnaseq_median_tumor_vaf_exp = median(readcounts[i,"RNAseq_Tumor_VAF"])
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor read coverage", rnaseq_median_tumor_coverage, "RNAseq","SNV", "Median", "Median tumor coverage (reads supporting reference and variant bases)")
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor read coverage - observed variants only", rnaseq_median_tumor_coverage_obs, "RNAseq","SNV", "Median", "Median tumor coverage (reads supporting reference and variant bases) of variants with at least one RNAseq read")  
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor read coverage - expressed genes only", rnaseq_median_tumor_coverage_exp, "RNAseq","SNV", "Median", "Median tumor coverage (reads supporting reference and variant bases) of variants in a gene that is expressed")
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor VAF", rnaseq_median_tumor_vaf, "RNAseq","SNV", "Median", "Median tumor variant allele frequency")
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor VAF - observed variants only", rnaseq_median_tumor_vaf_obs, "RNAseq","SNV", "Median", "Median tumor variant allele frequency of variants with at least one RNAseq read")  
  stats[dim(stats)[1]+1,] = c("RNAseq median tumor VAF - expressed genes only", rnaseq_median_tumor_vaf_exp, "RNAseq","SNV", "Median", "Median tumor variant allele frequency of variants in a gene that is expressed")
}
# - RNAseq Normal (all, variant observed only, expressed gene only)
if (length(which(names(readcounts)=="RNAseq_Normal_var_rc")) & length(which(names(readcounts)=="RNAseq_Normal_gene_FPKM"))){
  rnaseq_median_normal_coverage = median(readcounts[,"RNAseq_Normal_ref_rc"] + readcounts[,"RNAseq_Normal_var_rc"])
  rnaseq_median_normal_vaf = median(readcounts[,"RNAseq_Normal_VAF"])
  i = which(readcounts[,"RNAseq_Normal_var_rc"] > 0)
  rnaseq_median_normal_coverage_obs = median(readcounts[i,"RNAseq_Normal_ref_rc"] + readcounts[i,"RNAseq_Normal_var_rc"])
  rnaseq_median_normal_vaf_obs = median(readcounts[i,"RNAseq_Normal_VAF"])
  i = which(readcounts[,"RNAseq_Normal_gene_FPKM"] > fpkm_cutoff)
  rnaseq_median_normal_coverage_exp = median(readcounts[i,"RNAseq_Normal_ref_rc"] + readcounts[i,"RNAseq_Normal_var_rc"])
  rnaseq_median_normal_vaf_exp = median(readcounts[i,"RNAseq_Normal_VAF"])
  stats[dim(stats)[1]+1,] = c("RNAseq median normal read coverage", rnaseq_median_normal_coverage, "RNAseq","SNV", "Median", "Median normal coverage (reads supporting reference and variant bases)")
  stats[dim(stats)[1]+1,] = c("RNAseq median normal read coverage - observed variants only", rnaseq_median_normal_coverage_obs, "RNAseq","SNV", "Median", "Median normal coverage (reads supporting reference and variant bases) of variants with at least one RNAseq read")  
  stats[dim(stats)[1]+1,] = c("RNAseq median normal read coverage - expressed genes only", rnaseq_median_normal_coverage_exp, "RNAseq","SNV", "Median", "Median normal coverage (reads supporting reference and variant bases) of variants in a gene that is expressed")
  stats[dim(stats)[1]+1,] = c("RNAseq median normal VAF", rnaseq_median_normal_vaf, "RNAseq","SNV", "Median", "Median normal variant allele frequency")
  stats[dim(stats)[1]+1,] = c("RNAseq median normal VAF - observed variants only", rnaseq_median_normal_vaf_obs, "RNAseq","SNV", "Median", "Median normal variant allele frequency of variants with at least one RNAseq read")  
  stats[dim(stats)[1]+1,] = c("RNAseq median normal VAF - expressed genes only", rnaseq_median_normal_vaf_exp, "RNAseq","SNV", "Median", "Median normal variant allele frequency of variants in a gene that is expressed")
}

#Summarize variants in various categories:
#'supported' variants                -> (RNA-seq Variant Read Count >= 1)
#'missed' variants                   -> (RNA-seq Variant Read Count == 0) && (RNA-seq Total Read Count < rc_cutoff) && (RNA-seq Gene FPKM >= fpkm_cutoff)
#'expressed' variants.               -> (RNA-seq Total Read Count >= rc_cutoff) && (RNA-seq Variant Read Count >=1) && (RNA-seq Gene FPKM >= fpkm_cutoff)
#'mutant allele biased' variants.    -> ('expressed') && (VAF Difference > +vaf_diff_cutoff)
#'wild type allele biased' variants. -> ('expressed') && (VAF Difference < -vaf_diff_cutoff)
#'silent gene' variants              -> RNA-seq Read Count < rc_cutoff && Gene FPKM < fpkm_cutoff

#Description of categories
#'supported' variants                -> At least the minimal support for the variant found in RNA-seq
#'missed' variants                   -> Gene was expressed but variant position was not well covered AND mutation was not observed above minimal level
#'expressed' variants.               -> Variant had at least minimal support AND the variant position was well covered AND the gene was expressed
#'mutant allele biased' variants.    -> Variant position was well covered in WGS/Exome AND Variant position was well covered in RNA-seq AND Gene was expressed AND Variant allele frequency in RNA is skewed towards mutant allele
#'wild type allele biased' variants. -> Variant position was well covered in WGS/Exome AND Variant position was well covered in RNA-seq AND Gene was expressed AND Variant allele frequency in RNA is skewed towards wild type allele
#'silent gene' variants              -> Variant position was not well covered AND Gene was not expressed 

min_support_cutoff_text = paste(" (Variant Min. Support Cutoff = ", min_support_cutoff, ") ", sep="")
rc_cutoff_text = paste(" (Variant RC Cutoff = ", rc_cutoff, ") ", sep="")
fpkm_cutoff_text = paste(" (Gene FPKM Cutoff = ", fpkm_cutoff, ") ", sep="")
vaf_diff_cutoff_pos_text = paste(" (Variant allele frequency diff cutoff = ", vaf_diff_cutoff_pos, ") " , sep="")
vaf_diff_cutoff_neg_text = paste(" (Variant allele frequency diff cutoff = ", vaf_diff_cutoff_neg, ") " , sep="")

supported_variants_desc = paste("Variants from WGS/Exome at least minimal support in RNA-seq", min_support_cutoff_text, sep="")
expressed_variants_desc = paste("Variants from WGS/Exome exceeding a read cutoff with RNA-seq support", min_support_cutoff_text, rc_cutoff_text, fpkm_cutoff_text, sep="")
missed_variants_desc = paste("Variants from WGS/Exome NOT detected in RNAseq even though the gene was expressed", min_support_cutoff_text, rc_cutoff_text, fpkm_cutoff_text, sep="")
silent_gene_variants_desc = paste("Variants that appear to be in a gene that is not expressed", rc_cutoff_text, fpkm_cutoff_text, sep="")
mutant_biased_variants_desc = paste("Variants that appear to be preferentially expressed from the mutant allele", rc_cutoff_text, fpkm_cutoff_text, vaf_diff_cutoff_pos_text, sep="")
wildtype_biased_variants_desc = paste("Variants that appear to be preferentially expressed from the wild type allele", rc_cutoff_text, fpkm_cutoff_text, vaf_diff_cutoff_neg_text, sep="")

readcounts[,"MutantExpressionClass"] = "Other"
if (length(which(names(readcounts)=="RNAseq_Tumor_var_rc")) & length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  i = which(readcounts[,"RNAseq_Tumor_var_rc"] >= min_support_cutoff)
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "Supported"}
  supported_variants = length(i)
  supported_variants_p = round(((supported_variants/total_variants)*100), digits=2)

  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) >= rc_cutoff) & readcounts[,"RNAseq_Tumor_var_rc"] >= min_support_cutoff & readcounts[,"RNAseq_Tumor_gene_FPKM"] >= fpkm_cutoff)
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "Expressed"}
  expressed_variants = length(i)
  expressed_variants_p = round(((expressed_variants/total_variants)*100), digits=2)

  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) < rc_cutoff) & readcounts[,"RNAseq_Tumor_var_rc"] < min_support_cutoff & readcounts[,"RNAseq_Tumor_gene_FPKM"] >= fpkm_cutoff)
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "Missed"}
  missed_variants = length(i)
  missed_variants_p = round(((missed_variants/total_variants)*100), digits=2)

  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) < rc_cutoff) & readcounts[,"RNAseq_Tumor_gene_FPKM"] < fpkm_cutoff)
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "SilentGene"}
  silent_gene_variants = length(i)
  silent_gene_variants_p = round(((silent_gene_variants/total_variants)*100), digits=2)

  stats[dim(stats)[1]+1,] = c("Number of supported variants", supported_variants, "RNA-seq", "SNV", "Count", supported_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of supported variants", supported_variants_p, "RNA-seq", "SNV", "Percent", supported_variants_desc)

  stats[dim(stats)[1]+1,] = c("Number of expressed variants", expressed_variants, "RNA-seq", "SNV", "Count", expressed_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of expressed variants", expressed_variants_p, "RNA-seq", "SNV", "Percent", expressed_variants_desc)

  stats[dim(stats)[1]+1,] = c("Number of missed variants", missed_variants, "RNA-seq", "SNV", "Count", missed_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of missed variants", missed_variants_p, "RNA-seq", "SNV", "Percent", missed_variants_desc)

  stats[dim(stats)[1]+1,] = c("Number of variants in non-expressed genes", silent_gene_variants, "RNA-seq", "SNV", "Count", silent_gene_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of variants in non-expressed genes", silent_gene_variants_p, "RNA-seq", "SNV", "Percent", silent_gene_variants_desc)
}

#WGS vs RNAseq
if (length(which(names(readcounts)=="RNAseq_Tumor_VAF")) & length(which(names(readcounts)=="WGS_Tumor_VAF")) & length(which(names(readcounts)=="RNAseq_Tumor_gene_FPKM"))){
  #Calculate the variant allele frequency difference between WGS and RNAseq
  readcounts[,"WGS_Tumor_VAF_Diff"]=(readcounts[,"RNAseq_Tumor_VAF"] - readcounts[,"WGS_Tumor_VAF"])
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) >= rc_cutoff) & ((readcounts[,"WGS_Tumor_ref_rc"] + readcounts[,"WGS_Tumor_var_rc"]) >= rc_cutoff) & (readcounts[,"RNAseq_Tumor_gene_FPKM"] >= fpkm_cutoff) & (readcounts[,"WGS_Tumor_VAF_Diff"] > vaf_diff_cutoff_pos))
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "MutantBiased"}
  mutant_biased_variants = length(i)
  mutant_biased_variants_p = round(((mutant_biased_variants/total_variants)*100), digits=2)
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) >= rc_cutoff) & ((readcounts[,"WGS_Tumor_ref_rc"] + readcounts[,"WGS_Tumor_var_rc"]) >= rc_cutoff) & (readcounts[,"RNAseq_Tumor_gene_FPKM"] >= fpkm_cutoff) & (readcounts[,"WGS_Tumor_VAF_Diff"] < vaf_diff_cutoff_neg))
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "WildTypeBiased"}
  wildtype_biased_variants = length(i)
  wildtype_biased_variants_p = round(((wildtype_biased_variants/total_variants)*100), digits=2)

  stats[dim(stats)[1]+1,] = c("Number of WGS variants with mutant biased expression", mutant_biased_variants, "RNA-seq and WGS", "SNV", "Count", mutant_biased_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of WGS variants with mutant biased expression", mutant_biased_variants_p, "RNA-seq and WGS", "SNV", "Percent", mutant_biased_variants_desc)

  stats[dim(stats)[1]+1,] = c("Number of WGS variants with wild-type biased expression", wildtype_biased_variants, "RNA-seq and WGS", "SNV", "Count", wildtype_biased_variants_desc)    
  stats[dim(stats)[1]+1,] = c("Percent of WGS variants with wild-type biased expression", wildtype_biased_variants_p, "RNA-seq and WGS", "SNV", "Percent", wildtype_biased_variants_desc)    
}

#Exome vs RNAseq
if (length(which(names(readcounts)=="RNAseq_Tumor_VAF")) & length(which(names(readcounts)=="Exome_Tumor_VAF"))){
  #Calculate the variant allele frequency difference between WGS and RNAseq 
  readcounts[,"Exome_Tumor_VAF_Diff"]=(readcounts[,"RNAseq_Tumor_VAF"] - readcounts[,"Exome_Tumor_VAF"])
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) > rc_cutoff) & ((readcounts[,"Exome_Tumor_ref_rc"] + readcounts[,"Exome_Tumor_var_rc"]) > rc_cutoff) & (readcounts[,"RNAseq_Tumor_gene_FPKM"] > fpkm_cutoff) & (readcounts[,"Exome_Tumor_VAF_Diff"] > vaf_diff_cutoff_pos))
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "MutantBiased"}
  mutant_biased_variants = length(i)
  mutant_biased_variants_p = round(((mutant_biased_variants/total_variants)*100), digits=2)
  i = which(((readcounts[,"RNAseq_Tumor_ref_rc"] + readcounts[,"RNAseq_Tumor_var_rc"]) > rc_cutoff) & ((readcounts[,"Exome_Tumor_ref_rc"] + readcounts[,"Exome_Tumor_var_rc"]) > rc_cutoff) & (readcounts[,"RNAseq_Tumor_gene_FPKM"] > fpkm_cutoff) & (readcounts[,"Exome_Tumor_VAF_Diff"] < vaf_diff_cutoff_neg))
  if(length(i)){readcounts[i,"MutantExpressionClass"] = "WildTypeBiased"}
  wildtype_biased_variants = length(i)
  wildtype_biased_variants_p = round(((wildtype_biased_variants/total_variants)*100), digits=2)
  stats[dim(stats)[1]+1,] = c("Number of Exome variants with mutant biased expression", mutant_biased_variants, "RNA-seq and Exome", "SNV", "Count", mutant_biased_variants_desc)
  stats[dim(stats)[1]+1,] = c("Percent of Exome variants with mutant biased expression", mutant_biased_variants_p, "RNA-seq and Exome", "SNV", "Percent", mutant_biased_variants_desc)
  
  stats[dim(stats)[1]+1,] = c("Number of Exome variants with wild-type biased expression", wildtype_biased_variants, "RNA-seq and Exome", "SNV", "Count", wildtype_biased_variants_desc)	
  stats[dim(stats)[1]+1,] = c("Percent of Exome variants with wild-type biased expression", wildtype_biased_variants_p, "RNA-seq and Exome", "SNV", "Percent", wildtype_biased_variants_desc)    

}


#Write out the stats to a tsv file
filename = "Stats.tsv"
write.table(stats, file=filename, sep="\t", row.names=FALSE, quote=FALSE)

#Write out the amended data.frame so that the user can identify Mutant biased gene expression etc.
filename = "VariantExpressionSummary.tsv"
write.table(readcounts, file=filename, sep="\t", row.names=FALSE, quote=FALSE)


