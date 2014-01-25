#!/usr/bin/env Rscript
#Written by Malachi Griffith

args = (commandArgs(TRUE))
infile = args[1]
combined_vaf_cols_string = args[2]
normal_vaf_cols_string = args[3]
tumor_day0_vaf_cols_string = args[4]
tumor_day30_vaf_cols_string = args[5]
target_name = args[6]
outdir = args[7]

combined_vaf_cols = strsplit(combined_vaf_cols_string, " ")[[1]]
normal_vaf_cols = as.numeric(strsplit(normal_vaf_cols_string, " ")[[1]])
tumor_day0_vaf_cols = as.numeric(strsplit(tumor_day0_vaf_cols_string, " ")[[1]])
tumor_day30_vaf_cols = as.numeric(strsplit(tumor_day30_vaf_cols_string, " ")[[1]])


#Define input variables
#infile="/gscmnt/sata132/techd/mgriffit/aml_trios/H_KA-174556/H_KA-174556_final_filtered_coding_clean.tsv"
#combined_vaf_cols = c("normal_day0_VAF","tumor_day0_VAF","tumor_day30_VAF")
#normal_vaf_cols = c(34,37,40)
#tumor_day0_vaf_cols = c(43,46,49)
#tumor_day30_vaf_cols = c(52,55,58)
#target_name="AML_RMG"
#outdir="/gscmnt/sata132/techd/mgriffit/aml_trios/H_KA-174556/"

if (length(args) < 7){
  message_text1 = "Required arguments missing: ./SnvIndelReport.pm.R ..."
  stop(message_text1)
}

#Load in raw data
data = read.table(infile, header=TRUE, as.is = 1:16, sep="\t")

#Create a function to plot VAFs from replicates
plot_vafs = function(file, x_label, main_label, gene_i){  
  gene_count = length(gene_i)
  pdf(file)
  #Find the max value observed and use that to limit the y-axis
  y_max = max(data[gene_i,c(normal_vaf_cols,tumor_day0_vaf_cols,tumor_day30_vaf_cols)]) + 20
  round(y_max, digits = -1)  
  if (y_max > 100){
    y_max = 100
  }
  if (y_max < 25){
    y_max = 25
  }
  par(mar=c(6,6,4,2)+0.1)
  par(mgp=c(5,1,0))
  par(font.lab=2)

  #Adjust point sizes and gene labels if there are too many genes
  point_cex=1
  names_cex=1
  if (gene_count >= 25){
    fraction = gene_count/35
    point_cex = point_cex/(gene_count/35)
    names_cex = point_cex/(gene_count/45)
  }
  if (point_cex < 0.1){
    point_cex = 0.1
  }

  boxplot(t(data[gene_i,34]), names=data[gene_i,"default_gene_name"], ylim=c(0,y_max), border="white", color="white", xlab=x_label, ylab="VAF", main=main_label, las=2, xaxt="n")

  #Drop the gene labels entirely if there are too many genes
  if (gene_count < 50){
    axis(1, cex.axis=names_cex, labels=data[gene_i,"default_gene_name"], las=2, at=1:length(gene_i))
  }

  xp = 0
  for (i in gene_i){
    xp = xp+1
    flank = 0.125
   
    #Only display individual replicates if there not too many genes
    if (gene_count <= 25){
      xpoints = runif(n=length(normal_vaf_cols), min=(xp-flank), max=(xp+flank))
      points(x=xpoints, y=data[i,normal_vaf_cols], col="dark green", pch=16, cex=point_cex)
      points(x=xpoints, y=data[i,tumor_day0_vaf_cols], col="red", pch=16, cex=point_cex)
      points(x=xpoints, y=data[i,tumor_day30_vaf_cols], col="orange", pch=16, cex=point_cex)
      lines(x=c(xp-0.25,xp+0.25), y=rep(data[i,"normal_day0_VAF"],2), col="dark green", lwd=2)
      lines(x=c(xp-0.25,xp+0.25), y=rep(data[i,"tumor_day0_VAF"],2), col="red", lwd=2)
      lines(x=c(xp-0.25,xp+0.25), y=rep(data[i,"tumor_day30_VAF"],2), col="orange", lwd=2)
    }else{
      points(x=xp, y=data[i,"normal_day0_VAF"], col="dark green", pch=16, cex=point_cex)
      points(x=xp, y=data[i,"tumor_day0_VAF"], col="red", pch=16, cex=point_cex)
      points(x=xp, y=data[i,"tumor_day30_VAF"], col="orange", pch=16, cex=point_cex)
    }
  }
  legend("topright", legend=c("normal", "day0 tumor", "day30 tumor"), col=c("dark green","red","orange"), pch=16)
  dev.off()
}



#################################################################################################
#First plot the target genes only

#Find genes that overlapped the target gene list (e.g. RMG list)
target_gene_i = which(data[,target_name] == 1)
if (length(target_gene_i) > 0){
  outfile = paste(outdir, "target_gene_vafs.pdf", sep="")
  x_label = paste(target_name, " Gene Mutations", sep="")
  main_label = paste("Replicate VAFs for target genes: ", target_name, sep="")
  plot_vafs(outfile, x_label, main_label, target_gene_i)
}


#Now plot the top 10 VAF genes
o = order(data[,"max_tumor_vaf_observed"], decreasing=TRUE)
target_gene_i = o
if (length(o) >= 10){
  target_gene_i = o[1:10]
}

if (length(target_gene_i) > 0){
  outfile = paste(outdir, "top_gene_vafs.pdf", sep="")
  x_label = "Gene Mutations"
  main_label = "Replicate VAFs for highest VAF mutations"
  plot_vafs(outfile, x_label, main_label, target_gene_i)

}

#Now plot all genes
o = order(data[,"max_tumor_vaf_observed"], decreasing=TRUE)
target_gene_i = o
if (length(target_gene_i) > 0){
  outfile = paste(outdir, "all_gene_vafs.pdf", sep="")
  x_label = "Gene Mutations"
  main_label = "Replicate VAFs for all mutations"
  if (length(target_gene_i) > 25){
    main_label = "Aggregate VAFs for all mutations"
  }
  plot_vafs(outfile, x_label, main_label, target_gene_i)
}

