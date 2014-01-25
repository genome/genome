library(ggplot2)
library(plyr)

#This script takes a list of druggable genes broken down by patient and data type and plots visualization
druggable_genes_file="/Users/ogriffit/Dropbox/WashU/Projects/G4C/Clinseq_vs_druggable/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_collapsed.txt"
setwd("/Users/ogriffit/Dropbox/WashU/Projects/G4C/Clinseq_vs_druggable/")

data=read.table(file=druggable_genes_file, header=TRUE, sep="\t", as.is=c(1:6))
#Fix patient names
data[which(data[,"patient"]=="AML103"),"patient"]="ALL1"

interactions=data[,c(1,3,4)]
counts_table=table(interactions[,1:2])
counts=ddply(interactions, .(patient, data_type), function(d) length(unique(d$gene)))

#Change column and variable names for nice plotting
colnames(counts)=c("Patient","Type","Frequency")
counts[which(counts[,"Type"]=="exome"),"Type"]="Exome-Seq"
counts[which(counts[,"Type"]=="wgs"),"Type"]="WGS"
counts[which(counts[,"Type"]=="rnaseq"),"Type"]="RNA-Seq"

#Plot with actual counts
pdf(file="Clinseq_Druggable_Genes_by_Type_SM_DB_counts.pdf")
ggplot(counts,aes(x=factor(Patient),y=Frequency,fill=factor(Type))) + geom_bar(position="stack") + labs(x = "Patient", y = "Actionable Drug Targets") + labs(fill="") + theme_bw() + opts(title="Clin-Seq Predictions by Data Type", plot.title=theme_text(face="bold", size=20, hjust=1, vjust=1), axis.title.x=theme_text(face="bold", size=20), axis.title.y=theme_text(face="bold", size=20, angle=90, hjust=0.7), axis.text.x = theme_text(angle=90, hjust=1.2, size=16), axis.text.y = theme_text(angle=90, size=16, vjust=0.7), legend.text = theme_text(size=14))
dev.off()

#Plot with relative counts
pdf(file="Clinseq_Druggable_Genes_by_Type_SM_DB_fractions.pdf")
ggplot(counts,aes(x=factor(Patient),y=Frequency,fill=factor(Type))) + geom_bar(position="fill") + labs(x = "Patient", y = "Fraction Actionable Drug Targets") + labs(fill="") + theme_bw() + opts(title="Clin-Seq Predictions by Data Type", plot.title=theme_text(face="bold", size=20, hjust=1, vjust=1), axis.title.x=theme_text(face="bold", size=20), axis.title.y=theme_text(face="bold", size=20, angle=90, hjust=0.8), axis.text.x = theme_text(angle=90, hjust=1.2, size=16), axis.text.y = theme_text(angle=90, size=16, vjust=0.7), legend.text = theme_text(size=14))
dev.off()

#Plot pie chart of total unique number of genes per data type
exome_uniq_gene_count=length(unique(interactions[which(interactions[,"data_type"]=="exome"),"gene"]))
wgs_uniq_gene_count=length(unique(interactions[which(interactions[,"data_type"]=="wgs"),"gene"]))
rnaseq_uniq_gene_count=length(unique(interactions[which(interactions[,"data_type"]=="rnaseq"),"gene"]))

data=data.frame(DataType = c("Exome-Seq","RNA-Seq","WGS"), size = c(exome_uniq_gene_count, rnaseq_uniq_gene_count, wgs_uniq_gene_count))
data$p<-(p<-cumsum(data$size)-diff(c(0,cumsum(data$size)))*(1-0.5))
pdf(file="Clinseq_Druggable_Genes_by_Type_SM_DB_pie.pdf")
ggplot(data, aes(x = factor(1), fill = DataType, weight=size)) +  geom_bar(width = 1) + coord_polar(theta="y") + geom_text(aes(x= 1.6, y=p, angle=-p*360,label=size),vjust=0) + theme_bw() + opts(axis.title.x = theme_blank(), axis.title.y= theme_blank(), axis.text.x=theme_blank(), axis.text.y=theme_blank(), panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(), panel.border = theme_blank(), axis.ticks=theme_blank())
dev.off()

