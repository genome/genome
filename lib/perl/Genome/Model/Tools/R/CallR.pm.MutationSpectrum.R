

#This file contains code to graph a small region of interest from BamToCna data
plot_spectrum=function(spectrum_file="",output_file="",genome="",absolute_axis=T) {
    read.table(spectrum_file,fill=T,header=T,row.names=1)->spectrum;
    spectrum = (spectrum[0:6,]/spectrum['SNVs',1])*100;
    spectrum$Synonomous=c();
    spectrum=as.matrix(t(spectrum));
    pdf(file=output_file,width=6,height=6);
    title = "Mutation Spectrum";
    if(genome!="") {
        title = paste(title,"For",genome);
    }
    if(absolute_axis) {
        barplot(spectrum,beside=T,xlab="Mutation Class",ylim=c(0,100),ylab="Percent of Total Mutations",main=title,col=c("darkorange","darkgreen","purple4","darkred","darkblue","tan4"),space=c(0,0.1));
    }
    else {
        barplot(spectrum,beside=T,xlab="Mutation Class",ylab="Percent of Total Mutations",main=title,col=c("darkorange","darkgreen","purple4","darkred","darkblue","tan4"),space=c(0,0.1));
    }
    dev.off();
}

sequence_context_pvalue = function(p0=NULL,ts_before_cs=NULL,n=NULL,outfile=NULL,cs_before_cs=NULL,before_or_after="before") {
    p0 = as.numeric(p0);
    ts_before_cs = as.numeric(ts_before_cs);
    cs_before_cs = as.numeric(cs_before_cs);
    n = as.numeric(n);
    before_or_after = as.character(before_or_after);
    pyramidine_base_before_Cs_percentage = sprintf('%.2f',((ts_before_cs + cs_before_cs)*100)/n);
    # this is a proportion test. z-score followed by standard normal distribution p-value. #
    # sample_successes = ts_before_cs;
    sample_successes = ts_before_cs + cs_before_cs;
    pbar = sample_successes / n;
    z=(pbar-p0)/sqrt(p0*(1-p0)/n);
    pval = 2*pnorm(-abs(z),lower.tail=TRUE);

    # print results to file #
    cat(paste("Null hypothesis: Sample proportion is equivalent to hypothesized proportion.",sep=""),file=outfile,sep="\n");
    cat(paste("Hypothesized Porportion (Control Proportion): ",sprintf('%.3f',p0),sep=""),file=outfile,sep="\n",append=T);
    cat(paste("Sample Size: ",n,sep=""),file=outfile,sep="\n",append=T);
    cat(paste("Num of Ts ",before_or_after," Cs: ",ts_before_cs,sep=""),file=outfile,sep="\n",append=T);
    cat(paste("Num of Cs ",before_or_after," Cs: ",cs_before_cs,sep=""),file=outfile,sep="\n",append=T);
    cat(paste("% pyramidine bases ",before_or_after," Cs: ",pyramidine_base_before_Cs_percentage,"%",sep=""),file=outfile,sep="\n",append=T);
    # cat(paste("Sample Successes (Here, Ts before Cs): ",sample_successes,sep=""),file=outfile,sep="\n",append=T);
    cat(paste("Sample Successes (Here, Cs and Ts ",before_or_after," Cs): ",sample_successes,sep=""),file=outfile,sep="\n",append=T);
    cat(paste("Sample Proportion: ",sprintf('%.3f',pbar),sep=""),file=outfile,sep="\n",append=T);
    cat(paste("P-value: ",sprintf('%.5e',pval),sep=""),file=outfile,sep="\n",append=T);
}


plot_multi_mutation_spectrum <- function(input_file,plot_title="Mutation Spectrum",output_file=NULL,plot_type="facet1",pvalue=FALSE,plot.sample.order=NULL) {

  library("ggplot2");
  mutation_spectrum <- read.table(input_file,sep="\t",header=T);

  #mutation_spectrum$Sample <- toupper(mutation_spectrum$Sample); #convert sample label to uppercase
  mutation_spectrum$Density <- mutation_spectrum$Density/100;
  #separate out the basechange from transiton/transverstion
  mut_spec1 <- subset(mutation_spectrum,grepl("->",Category));
  mut_spec1$Category <- factor(mut_spec1$Category);
  mut_spec2 <- subset(mutation_spectrum,grepl("Tran",Category));
  mut_spec2$Category <- factor(mut_spec2$Category);

  if((pvalue)) {
    mut_spec1_pvalue <- calc_pvalue(mutation_spectrum,"->");
    mut_spec2_pvalue <- calc_pvalue(mutation_spectrum,"Tran");
  }else {
    mut_spec1_pvalue <- NULL;
    mut_spec2_pvalue <- NULL;
  }
  
  if(plot_type == "facet1") {
    p1 <- barplot_facet_mutation_type(data=mut_spec1,plot_title,pvalue=mut_spec1_pvalue,plot_order=plot.sample.order);
    p2 <- barplot_facet_mutation_type(data=mut_spec2,plot_title,pvalue=mut_spec2_pvalue,plot_order=plot.sample.order);
  }else if(plot_type == "facet2"){
    p1 <- barplot_facet_sample(data=mut_spec1,plot_title);
    p2 <- barplot_facet_sample(data=mut_spec2,plot_title);
  } else if(plot_type == "bar1") {
    p1 <- make_dodge_barplot(mut_spec1,plot_title,plot_order=plot.sample.order);
    p2 <- make_dodge_barplot(mut_spec2,plot_title,plot_order=plot.sample.order);
  } else if(plot_type == "bar2") {
    p1 <- make_stack_barplot(mut_spec1,plot_title,plot_order=plot.sample.order);
    p2 <- make_stack_barplot(mut_spec2,plot_title,plot_order=plot.sample.order);
  }
  
  #if outputfile is defined, print plots to outputfile
  #else, return ggplot_obj to caller
  if(!is.null(output_file)) {
    pdf(file=output_file,width=12,height=8);
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1,3)));
    print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=c(1,2)));
    print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=3));
    print(p1);
    print(p2);
    dev.off();
  }
  else {
    return(list(a=p1,b=p2));
  }

  
}


make_stack_barplot <- function(data.in,plot_title="Mutation Spectrum",plot_order=NULL) {

  #stack barplot with mutation category in different colors
  p <- ggplot(data.in,aes(x=Sample,y=Density));
  p <- p + geom_bar(position='fill',aes(fill=Category),width=0.9,stat="identity");
  p <- p + geom_bar(position='fill',aes(fill=Category),width=0.9, stat="identity", colour='gray22',legend=FALSE);
  if(length(levels(factor(data.in$Category))) > 2) {
    p <- p + scale_fill_brewer(pal="Set1",breaks=rev(levels(data.in$Category)));
  }else {
    p <- p + scale_fill_manual(value=c("red3","mediumblue"));
  }
 #allow user to specify the order of the X-axis
  if(!is.null(plot_order)) {
    x.axis.label.order <- strsplit(plot_order,',')[[1]];
    p <- p + scale_x_discrete(name="",limits=x.axis.label.order);
  }else {
    p <- p + scale_x_discrete(name="");
  }
  p <- p + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent");
  p <- p + opts(legend.position = 'right',legend.title=theme_blank());
  p <- p + opts(title=plot_title);
  plot_theme <- opts(axis.text.x=theme_text(colour='black',angle=90,hjust=1),axis.text.y=theme_text(colour='black'),plot.title=theme_text(size=14,face='bold'));
  p <- p + theme_bw() + plot_theme;

  return(p);

}

make_dodge_barplot <- function(data.in,plot_title="Mutation Spectrum",plot_order=NULL) {

  #change the order of the default factor of samples
  #this allows users to specify the order which each sample is plotted
  if(!is.null(plot_order)) {
    x.axis.label.order <- strsplit(plot_order,',')[[1]];
    data.in$Sample <- factor(data.in$Sample,levels=x.axis.label.order);   
  }

  p <- ggplot(data.in,aes(x=Category,y=Density));
  #dodge bar plot with samples in different colors
  p <- ggplot(data.in,aes(x=Category,y=Density));  
  p <- p + geom_bar(position='dodge',aes(fill=Sample),width=0.9,stat="identity");
  p <- p + geom_bar(position='dodge',aes(fill=Sample),width=0.9, stat="identity", colour='gray22',legend=FALSE);
  if(length(levels(factor(data.in$Sample))) > 2) {
    p <- p + scale_fill_hue();
  }else {
    p <- p + scale_fill_manual(value=c("red3","mediumblue"));
  }

  #specify the X and Y axis labels,legend,title,and plot theme
  p <- p + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent");
  #p <- p + scale_x_discrete(name="",breaks=c("A->C","A->G","A->T","C->A","C->G","C->T"), labels=c(" p = A1.5E-6","A->G0.002","A->T0.003","C->A0.004","C->G0.005","C->T0.006"));
  p <- p + scale_x_discrete(name="");
  p <- p + opts(legend.position = 'right',legend.title=theme_blank());
  p <- p + opts(title=plot_title);
  plot_theme <- opts(axis.text.x=theme_text(colour='black',angle=0,hjust=0.5),axis.text.y=theme_text(colour='black'),plot.title=theme_text(size=14,face='bold'));
  p <- p + theme_bw()+plot_theme;

  return(p);
  
}

barplot_facet_mutation_type <- function(data.in,plot_title="Mutation Spectrum",pvalue=NULL,plot_order=NULL) {

  #change the order of the default factor of samples
  #this allows users to specify the order which each sample is plotted
  if(!is.null(plot_order)) {
    x.axis.label.order <- strsplit(plot_order,',')[[1]];
    data.in$Sample <- factor(data.in$Sample,levels=x.axis.label.order);   
  }

  p <- ggplot(data.in,aes(y=Density,x=Category));
  p <- p + geom_bar(position='dodge',aes(fill=Sample),width=0.9);
  p <- p + geom_bar(position='dodge',aes(fill=Sample),width=0.9,colour='gray22',legend=FALSE);
  p <- p + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent");

  #pvalue =c("2.68e-11","< 2.2e-16","0.0401","< 2.2e-16","< 2.2e-16","< 2.2e-16");
  if(!is.null(pvalue)) {
    x.axis.breaks = names(pvalue);
    p <- p + scale_x_discrete(name="",breaks=x.axis.breaks, labels=pvalue);
  }else {
    p <- p + scale_x_discrete(name="");
  }
  p <- p + opts(legend.position = 'right',legend.title=theme_blank());

  if(length(levels(factor(data.in$Sample))) > 2) {
    p <- p + scale_fill_hue();
  }else {
    p <- p + scale_fill_manual(value=c("red3","mediumblue"));
  }
  p <- p + opts(title=plot_title);
  p <- p + facet_grid( . ~ Category,scales='free_x');

  plot_theme <- opts(axis.text.x=theme_text(colour='black'),axis.text.y=theme_text(colour='black'),plot.title=theme_text(size=14,face='bold'));
  p <- p + theme_bw()+ plot_theme;
  return(p);
}

barplot_facet_sample <- function(data,plot_title="Mutation Spectrum") {

  p <- ggplot(data,aes(y=Density,x=Category));
  p <- p + geom_bar(position='dodge',aes(fill=Category),width=0.9);
  p <- p + geom_bar(position='dodge',aes(fill=Category),width=0.9,colour='gray22',legend=FALSE);
  p <- p + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent");
  #p <- p + scale_x_discrete(name="",breaks=c("A->C","A->G","A->T","C->A","C->G","C->T"), labels=c(" p=1.5E-6","A->G0.002","A->T0.003","C->A0.004","C->G0.005","C->T0.006"));
  p <- p + scale_x_discrete(name="");
  p <- p + opts(legend.position = 'right',legend.title=theme_blank());

  if(length(levels(factor(data$Category))) > 2) {
    p <- p + scale_fill_hue();
  }else {
    p <- p + scale_fill_manual(value=c("red3","mediumblue"));
  }
  p <- p + opts(title=plot_title);
  p <- p + facet_wrap( ~ Sample,scales='free_x',nrow=1);

  plot_theme <- opts(axis.text.x=theme_text(colour='black',angle=90),axis.text.y=theme_text(colour='black'),plot.title=theme_text(size=14,face='bold'));
  p <- p + theme_bw()+plot_theme;
  
  return(p);
}


calc_pvalue <- function(inputDF,pattern="->") {

  #print(pattern);
  if(pattern=='Tran') {
    data2.in <- ddply(inputDF,.(Sample),subset,grepl("Tran",Category,fixed=TRUE));  #select only rows that matches pattern
  }else if(pattern=='->') {
    data2.in <- ddply(inputDF,.(Sample),subset,grepl("->",Category,fixed=TRUE)); 
  }else {
    data2.in <- ddply(inputDF,.(Sample),subset,grepl("->",Category,fixed=TRUE));  #select only rows that matches pattern
  }
  temp.df <- ddply(data2.in,.(Sample),numcolwise(sum)); #calc sum of mutations for each sample
  
  #make a named vector containing total number of mutations for each sample
  category.name <- as.character(temp.df$Sample);
  category.count <- temp.df$Count;
  names(category.count) <- category.name; 

  cat.label <- unique(as.character(data2.in$Category));
  pvalue_label=c();
  for (i in cat.label) {
    temp.df <- subset(data2.in,subset=(Category==i));
    temp.count <- temp.df$Count;
    names(temp.count) <- as.character(temp.df$Sample); 
    prop.test(temp.count,category.count) -> test.result;
    pvalue <- signif(test.result$p.value,digits=2);
    pvalue_label <- c(pvalue_label,pvalue); 
  }
  pvalue_label[which(pvalue_label<2.2e-16)]<-"p < 2.2e-16";
  pvalue_label[grep("p",pvalue_label,invert=T)] -> temp.pvalue;
  pvalue_label[grep("p",pvalue_label,invert=T)] <- paste("p =",temp.pvalue);
  names(pvalue_label) <- cat.label;
  #print(pvalue_label);
  
  return(pvalue_label);
}


make_dodge_barplot_facet_sample <- function (inputFile,plot_title="",num_row=NULL,outputFile=NULL,y_lim=c(0,1))  {
#takes the output of gmt analysis summarize mutation-spectrum
  
  library("ggplot2");
  
  mutation_spectrum <- read.table(inputFile,sep="\t",header=F);
  colnames(mutation_spectrum)=c("Category","count","percent","label")
  mutation_spectrum$label <- factor(mutation_spectrum$label);
  #separate out the basechange from transiton/transverstion
  mut_spec1 <- subset(mutation_spectrum,grepl("->",Category));
  mut_spec1$Category <- droplevels(mut_spec1$Category); #drop unused factor 
  mut_spec2 <- subset(mutation_spectrum,grepl("Tran",Category));
  mut_spec2$Category <- droplevels(mut_spec2$Category);
  
  p <- ggplot(mut_spec1,aes(x=Category,y=percent,fill=Category));
  p <- p + geom_bar(position='dodge',stat="identity",width=0.9);
  if (packageVersion("ggplot2") <= "0.8.9"){
    p <- p + geom_bar(position='dodge',stat="identity",width=0.9,colour='black',legend=FALSE);
  }else{
    p <- p + geom_bar(position='dodge',stat="identity",width=0.9,colour='black',show_guide=FALSE);
  }
  if(length(levels(mut_spec1$label)) > 1) { #turn on facet for multiple samples
    p <- p + facet_wrap( ~ label,nrow=num_row,scales='free_x')
  }
  p <- p + theme_bw();
  if (packageVersion("ggplot2") <= "0.8.9"){
    p <- p + scale_y_continuous(name='% Total Mutations', limits=c(0,1), formatter="percent");
  }else{
    library(scales)
    p <- p + scale_y_continuous(name='% Total Mutations', limits=c(0,1), labels=percent);
  }
  if (packageVersion("ggplot2") <= "0.8.9"){
    p <- p + scale_fill_brewer(pal="Dark2");
  }else{
    p <- p + scale_fill_brewer(palette="Dark2");
  }
  p <- p + scale_x_discrete(name="");
  p <- p + coord_cartesian(ylim=y_lim);
  if (packageVersion("ggplot2") <= "0.8.9"){
    p <- p + opts(title=plot_title,plot.title=theme_text(face="bold",size=16),legend.position = 'none',legend.title=theme_blank(),axis.text.x=theme_text(angle=90));
  }else{
    p <- p + labs(title=plot_title)
    p <- p + theme(plot.title=element_text(face="bold",size=16), legend.position = 'none', legend.title=element_blank(), axis.text.x=element_text(angle=90)); 
  }
  p2 <- ggplot(mut_spec2,aes(x=Category,y=percent,fill=Category));
  p2 <- p2 + geom_bar(position='dodge',stat="identity",width=0.9);
  if(length(levels(mut_spec2$label)) > 1) {  #turn on facet for multiple samples
    p2 <- p2 + facet_wrap( ~ label,nrow=num_row,scales='free_x')
  }
  p2 <- p2 + theme_bw();
  if (packageVersion("ggplot2") <= "0.8.9"){
    p2 <- p2 + scale_y_continuous(name='% Total Mutations', limits=c(0,1), formatter="percent");
  }else{
    library(scales)
    p2 <- p2 + scale_y_continuous(name='% Total Mutations', limits=c(0,1), labels=percent);
  }
  if (packageVersion("ggplot2") <= "0.8.9"){
    p2 <- p2 + scale_fill_manual(value=c("red3","mediumblue"));
  }else{
    p2 <- p2 + scale_fill_manual(values=c("red3","mediumblue"));
  }
  p2 <- p2 + scale_x_discrete(name="");
  if (packageVersion("ggplot2") <= "0.8.9"){ 
    p2 <- p2 + opts(title=plot_title,plot.title=theme_text(face="bold",size=16),legend.position = 'none',legend.title=theme_blank(),axis.text.x=theme_text(angle=90));
  }else{
    p2 <- p2 + labs(title=plot_title) 
    p2 <- p2 + theme(plot.title=element_text(face="bold",size=16), legend.position = 'none', legend.title=element_blank(), axis.text.x=element_text(angle=90)); 
  }
  if(!is.null(outputFile)) {
    pdf(file=outputFile,width=14,height=10);
    print(p);
    print(p2);
    dev.off();
  }
  else {
    return(list(a=p,b=p2));
  }
}


make_barplot_sample <- function (inputFile,plot_title="",label_angle=0,outputFile=NULL,y_lim=c(0,1),legend_bottom=FALSE,plot_type='dodge')  {

  library("ggplot2");
  
  mutation_spectrum <- read.table(inputFile,sep="\t",header=F);
  colnames(mutation_spectrum)=c("Category","count","percent","label")
  mutation_spectrum$label <- factor(mutation_spectrum$label);
  #separate out the basechange from transiton/transverstion
  mut_spec1 <- subset(mutation_spectrum,grepl("->",Category));
  #mut_spec1$Category <- factor(mut_spec1$Category);
  mut_spec1$Category <- droplevels(mut_spec1$Category); #drop unused factor 
  mut_spec2 <- subset(mutation_spectrum,grepl("Tran",Category));
  mut_spec2$Category <- droplevels(mut_spec2$Category);
  #mut_spec2$Category <- factor(mut_spec2$Category);

  p <- ggplot(mut_spec1,aes(x=label,y=percent,fill=Category));
  p <- p + geom_bar(position='dodge',stat="identity",width=0.9);
  p <- p + geom_bar(position='dodge',stat="identity",width=0.9,colour='black',legend=FALSE);
  if(plot_type == 'stack') { #switch to stack barplot
    p <- p + geom_bar(position='fill',stat="identity",width=0.9);
    p <- p + geom_bar(position='fill',stat="identity",width=0.9,colour='black',legend=FALSE);
    label_angle=90;  #if user specify stackbar graph, label should be turned 90 degrees
  }
  p <- p + theme_bw();
  p <- p + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent",expand=c(0,0));
  p <- p + scale_fill_brewer(pal="Dark2");
  p <- p + scale_x_discrete(name="",expand=c(0,0));
  #p <- p + coord_cartesian(ylim=y_lim);
  p <- p + opts(title=plot_title,plot.title=theme_text(face="bold",size=16),axis.text.x=theme_text(angle=label_angle));
  if(legend_bottom) {
    p <- p + opts(legend.position='bottom',legend.direction='horizontal',legend.title=theme_blank());
  }
  else {
    p <- p + opts(legend.position='right',legend.title=theme_blank());
  }
  
  p2 <- ggplot(mut_spec2,aes(x=label,y=percent,fill=Category));
  p2 <- p2 + geom_bar(position='dodge',stat="identity",width=0.9);
  if(plot_type == 'stack') { #switch to stack barplot
    p2 <- p2 + geom_bar(position='fill',stat="identity",width=0.9);
    p2 <- p2 + geom_bar(position='fill',stat="identity",width=0.9,colour='black',legend=FALSE);
    label_angle=90;  #if user specify stackbar graph, label should be turned 90 degrees 
  }
  p2 <- p2 + theme_bw();
  p2 <- p2 + scale_y_continuous(name='% Total Mutations',limits=c(0,1),formatter="percent");
  p2 <- p2 + scale_fill_manual(value=c("red3","mediumblue"));
  p2 <- p2 + scale_x_discrete(name="");
  p2 <- p2 + opts(title=plot_title,plot.title=theme_text(face="bold",size=16),axis.text.x=theme_text(angle=label_angle));
  if(legend_bottom) {
    p2 <- p2 + opts(legend.position='bottom',legend.direction='horizontal',legend.title=theme_blank());
  }
  else {
    p2 <- p2 + opts(legend.position='right',legend.title=theme_blank());
  }

  if(!is.null(outputFile)) {
    pdf(file=outputFile,width=14,height=10);
    print(p);
    print(p2);
    dev.off();
  }
  else {
    return(list(a=p,b=p2));
  }
  

}

plot_mutation_spectrum_seq_context <- function(input_file,plot_title=" ",output_file=NULL) {

  library("ggplot2");
  
  data.in=read.table(input_file,header=F,sep="\t");
  data.in$V1 = factor(data.in$V1,levels=c('A->C','C->A','A->G','C->G','A->T','C->T'));
  #plot_title="SMEL001003-D1-TB-01-0260";


  #stack barplot with mutation category in different colors
  p <- ggplot(data.in,aes(x=factor(V2),y=V4,fill=V3));
  p <- p + geom_bar(position='stack',stat="identity");
  p <- p + geom_bar(position='stack',stat="identity", colour='black',legend=FALSE);
  p <- p + facet_wrap( ~ V1, scales='free_y',nrow=3);
  #p <- p + scale_fill_brewer(pal="Set1",breaks=levels(data.in$V3));
  p <- p + scale_fill_manual(name="Base",value=c('A'=colors()[448],'T'=colors()[554],'G'=colors()[195],'C'=colors()[566]));
  p <- p + scale_x_discrete('Relative Position',expand=c(0,0));
  p <- p + scale_y_continuous("Number of Mutations",expand=c(0,0));
  p <- p + opts(title=plot_title,plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]),strip.text.x=theme_text(face="bold"),strip.background=theme_rect(fill=colors()[15]));

  if(!is.null(output_file)) {
    pdf(file=output_file,width=12,height=8);
    print(p);
    #print(p2);
    dev.off();
  }
  else {
    #return(list(a=p,b=p2));
    return(p);
  }
}

plot_mutation_spectrum_seq_contextV2 <- function(input4type,input2type,plot_title=" ",output_file=NULL) {
  #this version plots seq context for A/T/C/G and another plot for purine/pyrimidine

  library('ggplot2');
  
  data.in=read.table(input4type,header=F,sep="\t");
  data.in$V1 = factor(data.in$V1,levels=c('A','C','A->C','C->A','A->G','C->G','A->T','C->T'));

  #stack barplot with mutation category in different colors
  p1 <- ggplot(data.in,aes(x=factor(V2),y=V4,fill=V3));
  p1 <- p1 + geom_bar(position='stack',stat="identity");
  if (packageVersion("ggplot2") <= "0.8.9"){
    p1 <- p1 + geom_bar(position='stack',stat="identity", colour='black', legend=FALSE);
  }else{
    p1 <- p1 + geom_bar(position='stack',stat="identity", colour='black', show_guide=FALSE);
  }
  p1 <- p1 + facet_wrap( ~ V1, scales='free_y',nrow=4);
  if (packageVersion("ggplot2") <= "0.8.9"){
    p1 <- p1 + scale_fill_manual(name="Base",value=c('A'=colors()[448],'T'=colors()[554],'G'=colors()[195],'C'=colors()[566]));
  }else{
    p1 <- p1 + scale_fill_manual(name="Base",values=c('A'=colors()[448],'T'=colors()[554],'G'=colors()[195],'C'=colors()[566]));
  }
  p1 <- p1 + scale_x_discrete('Relative Position',expand=c(0,0));
  p1 <- p1 + scale_y_continuous("Number of Mutations",expand=c(0,0));
  if (packageVersion("ggplot2") <= "0.8.9"){
    p1 <- p1 + opts(title=plot_title,plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]),strip.text.x=theme_text(face="bold"),strip.background=theme_rect(fill=colors()[15]));
  }else{
    p1 <- p1 + labs(title=plot_title)
    p1 <- p1 + theme(plot.title = element_text(size=14, lineheight=.8, face="bold"), axis.text.y=element_text(colour='black'), axis.text.x=element_text(colour='black'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill=colors()[141]), strip.text.x=element_text(face="bold"), strip.background=element_rect(fill=colors()[15])); 
  }
  data.in2<-read.table(input2type,header=F,sep="\t");
  data.in2$V1 = factor(data.in2$V1,levels=c('A','C','A->C','C->A','A->G','C->G','A->T','C->T'));
  #stack barplot with mutation category in different colors
  p2 <- ggplot(data.in2,aes(x=factor(V2),y=V4,fill=V3));
  p2 <- p2 + geom_bar(position='stack',stat="identity");
  if (packageVersion("ggplot2") <= "0.8.9"){
    p2 <- p2 + geom_bar(position='stack',stat="identity", colour='black', legend=FALSE);
  }else{
    p2 <- p2 + geom_bar(position='stack',stat="identity", colour='black', show_guide=FALSE);
  }
  p2 <- p2 + facet_wrap( ~ V1, scales='free_y',nrow=4);
  if (packageVersion("ggplot2") <= "0.8.9"){
    p2 <- p2 + scale_fill_manual(name="Base",value=c('purine'=colors()[81],'pyrimidine'=colors()[93]));
  }else{
    p2 <- p2 + scale_fill_manual(name="Base",values=c('purine'=colors()[81],'pyrimidine'=colors()[93]));
  }
  p2 <- p2 + scale_x_discrete('Relative Position',expand=c(0,0));
  p2 <- p2 + scale_y_continuous("Number of Mutations",expand=c(0,0));
  if (packageVersion("ggplot2") <= "0.8.9"){
    p2 <- p2 + opts(title=plot_title,plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]),strip.text.x=theme_text(face="bold"),strip.background=theme_rect(fill=colors()[15]));
  }else{
    p2 <- p2 + labs(title=plot_title)
    p2 <- p2 + theme(plot.title = element_text(size=14, lineheight=.8, face="bold"), axis.text.y=element_text(colour='black'), axis.text.x=element_text(colour='black'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill=colors()[141]), strip.text.x=element_text(face="bold"),strip.background=element_rect(fill=colors()[15]));
  }
  #if output pdf file is not defined, return ggplot object
  #otherwise print the plot to the pdf file.
  if(!is.null(output_file)) {
    pdf(file=output_file,width=12,height=8);
    print(p1);
    print(p2);
    dev.off();
  }
  else {
    return(list(type4=p1,type2=p2));
  }
}

compare_prop2populations <- function(input_file,output_file) {

  input.df <- read.table(input_file,header=F,sep="\t");

  count_data <- input.df[,4:7];  #save the 4th-7th column into a separate dataframe
  p.values <- apply(count_data,1,do_proportion_test);
  p.values.adjusted <- p.adjust(p.values,method='fdr');

  input.df <- cbind(input.df,p.values,p.values.adjusted);
  
  write.table(input.df,file=output_file,quote=F,row.names=F,col.names=F,sep="\t");

}


do_proportion_test <- function(input_vector) {

  x <- input_vector[1:2];
  n <- input_vector[3:4];

  prop.test.obj <- prop.test(x,n);

  return(prop.test.obj$p.value);

}




plot_mutation_spectrum_bystrand <- function(inputFile,plot_title=NULL,outputFile=NULL,num_row=2,file_width=6,file_height=6 ) {

  library("ggplot2");

  data.in=read.table(inputFile,header=F,sep="\t");
  colnames(data.in) <- c("category","strand","count","label")
  data.in$category = factor(data.in$category,levels=c('A->C','C->A','A->G','C->G','A->T','C->T'));

  #plot_title="";

  p <- ggplot(data.in,aes(x=category,y=count,fill=strand));
  p <- p + geom_bar(position='dodge',stat="identity");
  p <- p + geom_bar(position='dodge',stat="identity", colour='black',legend=FALSE);
  if(length(levels(data.in$label)) > 1) { #turn on facet for multiple samples
    p <- p + facet_wrap( ~ label,nrow=num_row,scales='free')
  }
  p <- p + scale_fill_manual(name="Strand",value=c('transcribe'=colors()[448],'untranscribe'=colors()[554]),breaks=c('transcribe','untranscribe'),labels=c("Transcribed","Untranscribed"));
  p <- p + scale_x_discrete('');
  p <- p + scale_y_continuous("Number of Mutations");
  p <- p + opts(title=plot_title,plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.title.y=theme_text(angle=90,face='bold',size=14),axis.text.y=theme_text(colour='black',size=12),axis.text.x=theme_text(colour='black',size=12,angle=90),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]),strip.text.x=theme_text(face="bold",size=13),strip.background=theme_rect(fill=colors()[15]));
  
  if(!is.null(outputFile)) {
    pdf(file=outputFile,width=file_width,height=file_height);
    print(p);
    dev.off();
  }
  else {
    #return(list(a=p,b=p2));
    return(p);
  }



}


plot_mutation_distance<- function(mut_dist_file,chr_boundary_file,plot_title="Mutation Distance",output_file=NULL)  {

  library("ggplot2");
  mbp = 1e6;
  data.in <- read.table(mut_dist_file,sep="\t",header=F);
  data.in$V2 <- data.in$V2/mbp;
  colnames(data.in) <- c('chr','location','type','distance');

  #chr_boundary_file = "chr.length.merged.b37";
  chr_boundary <- read.table(chr_boundary_file,sep="\t",header=FALSE);
  colnames(chr_boundary)<- c('chr','start','end','xaxis_breakpt','label');  #add labels that make sense
  chr_boundary$start <- chr_boundary$start/mbp;
  chr_boundary$end <- chr_boundary$end/mbp;
  chr_boundary$xaxis_breakpt <-chr_boundary$xaxis_breakpt/1e6; 
  chr_boundary$chr[chr_boundary$chr==23] <- 'X'
  chr_boundary$chr[chr_boundary$chr==24] <- 'Y'
  
  p <- ggplot(data=chr_boundary);
  p <- p + geom_rect(aes(xmin=start,xmax=end,ymin=0.7,ymax=100000000,fill=label));
  p <- p + geom_point(data=data.in,aes(x=location,y=distance,colour=type),size=1.5,alpha=1);
  p <- p + geom_hline(data=data.in,yintercept=10000,colour=colors()[190],linetype='dashed');
  p <- p + geom_vline(aes(xintercept=end),colour='black');
  p <- p + scale_x_continuous(name='chromosome',expand=c(0,0),breaks=chr_boundary$xaxis_breakpt,labels=chr_boundary$chr);
  p <- p + scale_y_continuous(name='Intermutation Distance(bp)',trans='log10');
  p <- p <- p + coord_cartesian(ylim=c(0.7,1e8));
  p <- p + scale_colour_manual(name='',values=c('A->G'='orange','C->T'=colors()[97],'A->C'='red','C->G'=colors()[68],'C->A'=colors()[565],'A->T'=colors()[81]));
  p <- p + scale_fill_manual(name='',values=c('odd'='white','even'=colors()[235]),legend=FALSE);
  p <- p + opts(title=plot_title,plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.title.y=theme_text(angle=90,face='bold',size=11),axis.text.y=theme_text(colour='black',size=12),axis.title.x=theme_text(angle=0,face='bold',size=11),axis.text.x=theme_text(colour='black',size=8,angle=0),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),strip.text.x=theme_text(face="bold",size=10),panel.background=theme_rect(fill=colors()[2]),panel.border=theme_rect(colour='black')); 

  if(!is.null(output_file)) {
    pdf(file=output_file,width=14,height=5);
    #grid.newpage()
    #pushViewport(viewport(layout=grid.layout(1,3)));
    #print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=c(1,2)));
    #print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=3));
    print(p);
    #print(p2);
    dev.off();
  }
  else {
    #return(list(a=p1,b=p2));
    return(p);
  }

}
