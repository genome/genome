clonal_evolution = function(sample1,sample2,intersection_datafile,sample1_only_datafile,sample2_only_datafile,output_pdf,cluster=0)
{

#perform file tests
int_file_info = file.info(intersection_datafile);
sample1_file_info = file.info(sample1_only_datafile);
sample2_file_info = file.info(sample2_only_datafile);

#read in files
initialization_data = data.frame(V1=NA,V2=NA,V3=NA,V4=NA,V5=NA,V6=NA,V7=NA,V8=NA,V9=NA,V10=NA,V11=NA,V12=NA,V13=NA,V14=NA,V15=NA,V16=NA,V17=NA,V18=NA,V19=NA);
int = initialization_data;
sample1_only = initialization_data;
sample2_only = initialization_data;
if (int_file_info$size) { int = read.table(file=intersection_datafile,skip=0,sep="\t"); };
if (sample1_file_info$size) { sample1_only = read_varscan_file(sample1_only_datafile); };
if (sample2_file_info$size) { sample2_only = read_varscan_file(sample2_only_datafile); };

#create frequency variables
sample1_int_freqs = int$V3;
sample2_int_freqs = int$V4;
sample1_only_freqs = sample1_only$V11;
sample2_only_freqs = sample2_only$V11;
all_sample1_freqs = c(sample1_int_freqs,sample1_only_freqs);
all_sample2_freqs = c(sample2_int_freqs,sample2_only_freqs);

#create density plot variables
from = 0;
to = 100;
all_sample1_dens=NULL;
all_sample2_dens=NULL;
sample1_int_dens=NULL;
sample2_int_dens=NULL;
sample1_only_dens=NULL;
sample2_only_dens=NULL;

if (length(all_sample1_freqs)>2) { all_sample1_dens = density(all_sample1_freqs,from=from,to=to,na.rm=TRUE); }
if (length(all_sample2_freqs)>2) { all_sample2_dens = density(all_sample2_freqs,from=from,to=to,na.rm=TRUE); }
if (length(sample1_int_freqs)>1) { sample1_int_dens = density(sample1_int_freqs,from=from,to=to,na.rm=TRUE); }
if (length(sample2_int_freqs)>1) { sample2_int_dens = density(sample2_int_freqs,from=from,to=to,na.rm=TRUE); }
if (length(sample1_only_freqs)>1) { sample1_only_dens = density(sample1_only_freqs,from=from,to=to,na.rm=TRUE); }
if (length(sample2_only_freqs)>1) { sample2_only_dens = density(sample2_only_freqs,from=from,to=to,na.rm=TRUE); }

#plot these density diagrams
sample1 = as.character(gsub('_',' ',sample1));
sample2 = as.character(gsub('_',' ',sample2));
pdf(file=output_pdf,width=10,height=10);
par(mfrow=c(3,2));
if (!is.null(all_sample1_dens)) { plot_density(all_sample1_dens,sample1,"All SNVs");} else { plot_no_data(sample1,"All SNVs"); }
if (!is.null(all_sample2_dens)) { plot_density(all_sample2_dens,sample2,"All SNVs");} else { plot_no_data(sample2,"All SNVs"); }
if (!is.null(sample1_int_dens)) { plot_density(sample1_int_dens,sample1,"Shared SNVs");} else { plot_no_data(sample1,"Shared SNVs"); }
if (!is.null(sample2_int_dens)) { plot_density(sample2_int_dens,sample2,"Shared SNVs");} else { plot_no_data(sample2,"Shared SNVs"); }
if (!is.null(sample1_only_dens)) { plot_density(sample1_only_dens,sample1,paste(sample1,"-Only SNVs",sep=""));} else { plot_no_data(sample1,paste(sample1,"-Only SNVs",sep="")); }
if (!is.null(sample2_only_dens)) { plot_density(sample2_only_dens,sample2,paste(sample2,"-Only SNVs",sep=""));} else { plot_no_data(sample2,paste(sample2,"-Only SNVs",sep="")); }
dev.off();

}

read_varscan_file = function(filename,header=F,sep="\t")
{
    skip = NULL; if (header) { skip = 1; } else { skip = 0; }
    x<-read.table(file=filename,skip=skip,sep=sep);
    x$V7 = as.numeric(sub('%','',(x$V7)));
    x$V11 = as.numeric(sub('%','',(x$V11)));
    invisible(x);
}

plot_density = function(density,sample,title_string,color="black")
{
    plot(density$x,density$y,type="l",xlab=paste(sample," Variant Allele Frequency",sep=""),ylab=paste("Density of Mutations (N=",density$n,")",sep=""),col=color);
    title(title_string);
}

plot_no_data = function(sample,title_string,color="black")
{
    plot(0,0,xlab=paste(sample," Variant Allele Frequency",sep=""),ylab=paste("Density of Mutations (N=0)",sep=""),col=color);
    text(0,0,"NO DATA\n\n(N=0)",cex=2.0);
    title(title_string);
}
