####### GC Bias Library ###################################

gc_cov_bins = function (in.file,out.file,mean_read_length,mean_coverage=NULL) {

#read the input data
x = read.table(in.file,header=TRUE,sep="\t",quote="");
mean_read_length = as.numeric(mean_read_length);

#calculate the mean_coverage if not given
#Note: 1000bp windows are currently hard-coded in gc-bias tool bam-window call
if (is.null(mean_coverage)) { mean_coverage = mean(x$coverage) * mean_read_length / 1000; }
else { mean_coverage = as.numeric(mean_coverage); }

#calculate total reads
total_reads = sum(x$coverage);

#calculate per-bin coverage
bincov = c(mean(x$coverage[(x$gc_content>0)&(x$gc_content<=10)]),mean(x$coverage[(x$gc_content>10)&(x$gc_content<=20)]),mean(x$coverage[(x$gc_content>20)&(x$gc_content<=30)]),mean(x$coverage[(x$gc_content>30)&(x$gc_content<=40)]),mean(x$coverage[(x$gc_content>40)&(x$gc_content<=50)]),mean(x$coverage[(x$gc_content>50)&(x$gc_content<=60)]),mean(x$coverage[(x$gc_content>60)&(x$gc_content<=70)]),mean(x$coverage[(x$gc_content>70)&(x$gc_content<=80)]),mean(x$coverage[(x$gc_content>80)&(x$gc_content<=90)]),mean(x$coverage[(x$gc_content>90)&(x$gc_content<=100)]));

bincov = bincov * mean_read_length / 1000;

#calculate per-bin read-count
binreads = c(sum(x$coverage[(x$gc_content>0)&(x$gc_content<=10)]),sum(x$coverage[(x$gc_content>10)&(x$gc_content<=20)]),sum(x$coverage[(x$gc_content>20)&(x$gc_content<=30)]),sum(x$coverage[(x$gc_content>30)&(x$gc_content<=40)]),sum(x$coverage[(x$gc_content>40)&(x$gc_content<=50)]),sum(x$coverage[(x$gc_content>50)&(x$gc_content<=60)]),sum(x$coverage[(x$gc_content>60)&(x$gc_content<=70)]),sum(x$coverage[(x$gc_content>70)&(x$gc_content<=80)]),sum(x$coverage[(x$gc_content>80)&(x$gc_content<=90)]),sum(x$coverage[(x$gc_content>90)&(x$gc_content<=100)]));

#calculate normalized per-bin coverage
normbincov = bincov / mean_coverage;

#calculate normalized per-bin read-count
normbinreads = binreads / total_reads;

#remove 'NaN' values from coverage data - replace with zeros
bincov[is.na(bincov)]=0;
normbincov[is.na(normbincov)]=0;

#print some statistics for the data in the summary
cat("mean read-length: ",mean_read_length,"\n",file=out.file,fill=FALSE,append=TRUE);
cat("mean coverage: ",mean_coverage,"\n",file=out.file,fill=FALSE,append=TRUE);
cat("total reads: ",total_reads,"\n",file=out.file,fill=FALSE,append=TRUE);

#print un-normalized bin data
cat("\nNon-normalized coverage data per bin:\n",file=out.file,fill=FALSE,append=TRUE);
cat(paste("bin0-10","bin11-20","bin21-30","bin31-40","bin41-50","bin51-60","bin61-70","bin71-80","bin81-90","bin91-100\n",sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat(paste(bincov[1],bincov[2],bincov[3],bincov[4],bincov[5],bincov[6],bincov[7],bincov[8],bincov[9],bincov[10],sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat("\n",file=out.file,fill=FALSE,append=TRUE);

cat("\nNon-normalized read-count data per bin:\n",file=out.file,fill=FALSE,append=TRUE);
cat(paste("bin0-10","bin11-20","bin21-30","bin31-40","bin41-50","bin51-60","bin61-70","bin71-80","bin81-90","bin91-100\n",sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat(paste(binreads[1],binreads[2],binreads[3],binreads[4],binreads[5],binreads[6],binreads[7],binreads[8],binreads[9],binreads[10],sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat("\n",file=out.file,fill=FALSE,append=TRUE);

#print normalized bin data
cat("\nNormalized coverage data per bin:\n",file=out.file,fill=FALSE,append=TRUE);
cat(paste("bin0-10","bin11-20","bin21-30","bin31-40","bin41-50","bin51-60","bin61-70","bin71-80","bin81-90","bin91-100\n",sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat(paste(normbincov[1],normbincov[2],normbincov[3],normbincov[4],normbincov[5],normbincov[6],normbincov[7],normbincov[8],normbincov[9],normbincov[10],sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat("\n",file=out.file,fill=FALSE,append=TRUE);

cat("\nNormalized read-count data per bin:\n",file=out.file,fill=FALSE,append=TRUE);
cat(paste("bin0-10","bin11-20","bin21-30","bin31-40","bin41-50","bin51-60","bin61-70","bin71-80","bin81-90","bin91-100\n",sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat(paste(normbinreads[1],normbinreads[2],normbinreads[3],normbinreads[4],normbinreads[5],normbinreads[6],normbinreads[7],normbinreads[8],normbinreads[9],normbinreads[10],sep="\t"),file=out.file,fill=FALSE,append=TRUE);
cat("\n",file=out.file,fill=FALSE,append=TRUE);

}
