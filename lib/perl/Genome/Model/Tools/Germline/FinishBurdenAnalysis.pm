package Genome::Model::Tools::Germline::FinishBurdenAnalysis;

use warnings;
use strict;
use Carp;
use Genome;
use IO::File;
use POSIX qw( WIFEXITED );

class Genome::Model::Tools::Germline::FinishBurdenAnalysis {
  is => 'Genome::Model::Tools::Music::Base',
  has_input => [
    input_directory => { is => 'Text', doc => "Directory of Results of the Burden Analysis" },
    output_file => { is => 'Text', doc => "File with the Results of the Burden Analysis" },
    project_name => { is => 'Text', doc => "The name of the project", default => "Burden Analysis Results"},
    base_R_commands => { is => 'Text', doc => "The base R command library", default => '/gscuser/qzhang/gstat/burdentest/burdentest.R' },
  ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Combine results from a burden analysis on germline (PhenotypeCorrelation) data"                 
}

sub help_synopsis {
    return <<EOS
Run a burden analysis on germline (PhenotypeCorrelation) data
EXAMPLE:	gmt germline finish-burden-analysis --help
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
  return (
<<"EOS"
Run a burden analysis on germline (PhenotypeCorrelation) data
EXAMPLE:	gmt germline finish-burden-analysis --help
EOS
    );
}


###############

sub execute {                               # replace with real execution logic.
	my $self = shift;

    my $input_directory = $self->input_directory;
    my $base_R_commands = $self->base_R_commands;
    my $project_name = $self->project_name;

    my $output_file = $self->output_file;
    my $fh_outfile = new IO::File $output_file,"w";
    unless ($fh_outfile) {
        die "Failed to create output file $output_file!: $!";
    }

    opendir DIR, $input_directory or die "cannot open dir $input_directory: $!";
    my @files = grep { $_ ne '.' && $_ ne '..'} readdir DIR;
    closedir DIR;

    my @results_files;
    my @null_files;
    my @rarevariant_files;
    foreach my $file (@files) {
        if ($file =~ m/burden/) {
            push(@results_files,$file);
        }
        elsif ($file =~ m/null/) {
            push(@null_files,$file);
        }
        elsif ($file =~ m/single/) {
            push(@rarevariant_files,$file);
        }
        elsif ($file =~ m/err/) {
            die "Error file found $file, please resolve and/or re-run burden analysis for this variant-phenotype combination\n";
        }
    }


    foreach my $file (@null_files) {
        $file =~ s/.null$//;
        my ($pheno, $gene)  = split(/_/, $file);
        print "Null: $pheno\t$gene\n";
    }

    my %gene_variant_hash;
    foreach my $file (@rarevariant_files) {
        my $infh = new IO::File "$input_directory/$file","r";
        $file =~ s/.single.csv$//;
        my ($pheno, $gene)  = split(/_/, $file);
        my $header = $infh->getline;
        chomp($header);
        while (my $line = $infh->getline) {
            chomp($line);
            #Variant,Beta,SE,t,P
            my ($variant_name, $beta, $SE, $t, $P) = split(/,/, $line);
            my ($chr, $pos, $ref, $var) = split(/_/, $variant_name);
            unless ($chr =~ m/X/i || $chr =~ m/Y/i) {
                $chr =~ s/^\D+//;
            }
            else {
                $chr =~ s/^\D+X//;
                $chr =~ s/^\D+Y//;
            }
            my $chrpos = "$chr:$pos";
            $gene_variant_hash{$gene}{$chrpos}++;
        }
    }

    my $first = 1;
    foreach my $file (@results_files) {
        my $infh = new IO::File "$input_directory/$file","r";
        $file =~ s/.burden.csv$//;
        my ($pheno, $gene)  = split(/_/, $file);
        my $header = $infh->getline;
        chomp($header);
        if ($first) {
            print $fh_outfile "$header,Chromosome,Average_Position\n";
            $first = 0;
        }
        while (my $line = $infh->getline) {
            chomp($line);
            #heightRES,PSRC1,4384,5,0.01,0.958703922888546,0.9575,0.8092,0.9572,0.5614,0.6765,0.4057,0.4339
            my ($Trait,$Gene,$N,$V,$MAF,$CMC,$pCMC,$WSS,$aSum,$PWST,$SPWST,$SPWST_up,$SPWST_down) = split(/,/, $line);
            my @chrpositions = (sort keys %{$gene_variant_hash{$Gene}});
            my $average_position = 0;
            my $chr;
            foreach my $chrpos (@chrpositions) {
                ($chr, my $pos) = split(/:/, $chrpos);
                $average_position += $pos;
            }
            $average_position /= scalar(@chrpositions);
            print $fh_outfile "$line,$chr,$average_position\n";
        }
    }

    my $final_file = $output_file."_FDR";
    my $plot_file = $output_file.".pdf";

    my $R_burden_finisher_file = "$input_directory/Burden_Finisher.R";
    my $fh_R_finisher = new IO::File $R_burden_finisher_file,"w";
    unless ($fh_R_finisher) {
        die "Failed to create R options file $R_burden_finisher_file!: $!";
    }
    #-------------------------------------------------
    my $R_command_finisher = <<"_END_OF_R_";
missing.data=c("NA",".","");

x<-read.table("$output_file", sep = ",", header = TRUE);

x[x == 0] = 0.0000000001;

x\$CMC_log10=-log10(x\$CMC);
x\$pCMC_log10=-log10(x\$pCMC);
x\$WSS_log10=-log10(x\$WSS);
x\$aSum_log10=-log10(x\$aSum);
x\$PWST_log10=-log10(x\$PWST);
x\$SPWST_log10=-log10(x\$SPWST);
x\$SPWST.up_log10=-log10(x\$SPWST.up);
x\$SPWST.down_log10=-log10(x\$SPWST.down);
rownames(x) <- paste(x\$Trait,\"_\",x\$Gene,sep = \"\");

chromosomes <- sort(unique(x\$Chromosome));

#outfile_basename <- \"$output_file\";

#BEGIN PLOTTING IMAGE
pdf(file=\"$plot_file\",width=10,height=7.5,bg=\"white\");
colors <- rainbow(8);
library(genefilter);
for(k in 1:length(chromosomes)) {
    current_chr <- chromosomes[k];
#    outfile <- paste(outfile_basename,\"_chr\",current_chr,\".pdf\", sep = \"\");
    x_chr <- subset(x,x\$Chromosome == current_chr);
    avg_pos <- x_chr\$Average_Position;
    ymax <- max(x\$CMC_log10,x\$pCMC_log10,x\$WSS_log10,x\$aSum_log10,x\$PWST_log10,x\$SPWST_log10,x\$SPWST.up_log10,x\$SPWST.down_log10,na.rm = TRUE);

#subset chromosomes with wide regions
    max_avg_pos <- max(avg_pos);
    min_avg_pos <- min(avg_pos);
    diff_avg <- diff(sort(avg_pos));
    max_diff <- max(diff_avg);
    diff_position <- which(diff_avg==max_diff);
    sort_avg_pos <- sort(avg_pos)[diff_position];

    #library(genefilter)


    if (length(unique(avg_pos)) > 1) {
        if (max_diff > (0.5 * (max_avg_pos - min_avg_pos))){
            x_chr\$avg_pos_scaled <- genescale(avg_pos, axis=1, method=c(\"R\"), na.rm=TRUE)
print(x_chr\$avg_pos_scaled)[1:10];
            x_chr\$avg_pos_scaled <- x_chr[order(x_chr\$avg_pos_scaled),]
print(x_chr\$avg_pos_scaled)[1:10];
#            avg_pos_scaled <- genescale(avg_pos, axis=1, method=c(\"Z\"), na.rm=TRUE)
#            x_chr\$Average_Position <- avg_pos_scaled;
        }
    }



if(current_chr == 1) {
#print(x_chr\$Average_Position);
}
    gene_names <- unique((x_chr\$Gene));
    trait_names <- unique((x_chr\$Trait));

    #BEGIN PLOTTING -LOG10 P-VALUES
    for(j in 1:length(trait_names)) {
        current_trait <- trait_names[j];
        x_trait <- subset(x_chr,x_chr\$Trait == current_trait);
        avg_pos_trait <- x_trait\$Average_Position;
        avg_pos_range <- paste(round(min_avg_pos),\"-\",round(max_avg_pos),sep = \"\");
        plot(avg_pos_trait,x_trait\$CMC_log10,type="p",pch=16,cex=0.8,col=colors[1],main=paste(\"Burden Test Results for Chromosome \",current_chr,\":\",avg_pos_range,\" by Phenotype: \",current_trait,sep = \"\"),xlab=\"Average Chromosomal Position For the Rare SNVs (MAF<1%) in the Gene\",ylab=\"-log10\(p-value\)\",ylim=c(0,ymax*1.25), xaxt="n");
        points(avg_pos_trait,x_trait\$pCMC_log10,type="p",pch=16,cex=0.8,col=colors[2]);
        points(avg_pos_trait,x_trait\$WSS_log10,type="p",pch=16,cex=0.8,col=colors[3]);
        points(avg_pos_trait,x_trait\$aSum_log10,type="p",pch=16,cex=0.8,col=colors[4]);
        points(avg_pos_trait,x_trait\$PWST_log10,type="p",pch=16,cex=0.8,col=colors[5]);
        points(avg_pos_trait,x_trait\$SPWST_log10,type="p",pch=16,cex=0.8,col=colors[6]);
        points(avg_pos_trait,x_trait\$SPWST.up_log10,type="p",pch=16,cex=0.8,col=colors[7]);
        points(avg_pos_trait,x_trait\$SPWST.down_log10,type="p",pch=16,cex=0.8,col=colors[8]);

        plot_types <- c("CMC","pCMC","WSS","aSum","PWST","SPWST","SPWST.up","SPWST.down");
        legend(x=\"topright\", title = \"Burden Tests\", legend=plot_types,col=colors,pch=16,ncol=2);

        #Select -log10 > 1.5 for text labels below
        x2_trait <- subset(x_trait,x_trait\$CMC_log10 > 1.5 | x_trait\$pCMC_log10 > 1.5 | x_trait\$WSS_log10 > 1.5 | x_trait\$aSum_log10 > 1.5 | x_trait\$PWST_log10 > 1.5 | x_trait\$SPWST_log10 > 1.5 | x_trait\$SPWST.up_log10 > 1.5 | x_trait\$SPWST.down_log10 > 1.5, select=c(CMC_log10,pCMC_log10,WSS_log10,aSum_log10,PWST_log10,SPWST_log10,SPWST.up_log10,SPWST.down_log10));
        x3_trait <- subset(x_trait,x_trait\$CMC_log10 > 1.5 | x_trait\$pCMC_log10 > 1.5 | x_trait\$WSS_log10 > 1.5 | x_trait\$aSum_log10 > 1.5 | x_trait\$PWST_log10 > 1.5 | x_trait\$SPWST_log10 > 1.5 | x_trait\$SPWST.up_log10 > 1.5 | x_trait\$SPWST.down_log10 > 1.5, select=c(Average_Position));
        x4_trait <- subset(x_trait,x_trait\$CMC_log10 > 1.5 | x_trait\$pCMC_log10 > 1.5 | x_trait\$WSS_log10 > 1.5 | x_trait\$aSum_log10 > 1.5 | x_trait\$PWST_log10 > 1.5 | x_trait\$SPWST_log10 > 1.5 | x_trait\$SPWST.up_log10 > 1.5 | x_trait\$SPWST.down_log10 > 1.5, select=c(Gene));
        if(dim(x2_trait)[1]) {
            for(i in 1:length(x2_trait)) {
                if (!is.na(x2_trait[i]) && (x2_trait[i] > 2.5)) {
                    text(cbind(x3_trait,x2_trait[i]),lab=round(x2_trait[i],4), cex=0.7, pos=4);
                }
            }
        }

        # Make x axis tick marks without labels
        axis(1, at=unique(x_trait\$Average_Position), lab=F, tck=-0.008);
        # Plot x axis labels at default tick marks with labels at 45 degree angle
        text(unique(x_trait\$Average_Position), par("usr")[3], srt=45, adj = c(1.1,1.2), labels=unique(x_trait\$Gene), xpd=T, cex=1);

    }

    x_log10 <- subset(x, select=c(CMC_log10,pCMC_log10,WSS_log10,aSum_log10,PWST_log10,SPWST_log10,SPWST.up_log10,SPWST.down_log10));
}

#Q-Q Plots
for(i in 1:length(plot_types)) {
    index <- seq(1, nrow(x_log10[i]));
    uni <- index/nrow(x_log10[i]);
    loguni <- -log10(uni);
    x_plot <- t(sort(loguni));
    x_plot_label = paste(colnames(x_log10[i]),"_uniform", sep = \"\");

    y_plot <- sort(t(x_log10[i]));
    y_plot_label <- colnames(x_log10[i]);

    maxplot <- max(x_plot, y_plot, na.rm = TRUE);
    qqplot(x=x_plot,y=y_plot, xlab=x_plot_label,ylab=y_plot_label, xlim=c(0,maxplot),ylim=c(0,maxplot));
    abline(a=0,b=1);
}

for(i in 1:(length(plot_types) - 1)) {
    x_plot = x_log10[i];
    x_plot_label = colnames(x_log10[i]);
    for(j in (i+1):(length(plot_types))) {
        y_plot = x_log10[j];
        y_plot_label = colnames(x_log10[j]);
        maxplot <- max(x_plot, y_plot, na.rm = TRUE);
        qqplot(x=t(x_plot),y=t(y_plot), xlab=x_plot_label,ylab=y_plot_label, xlim=c(0,maxplot), ylim=c(0,maxplot));
        abline(a=0,b=1);
    }
}

dev.off();

x\$CMC_FDR=p.adjust(x\$CMC,method="fdr");
x\$pCMC_FDR=p.adjust(x\$pCMC,method="fdr");
x\$WSS_FDR=p.adjust(x\$WSS,method="fdr");
x\$aSum_FDR=p.adjust(x\$aSum,method="fdr");
x\$PWST_FDR=p.adjust(x\$PWST,method="fdr");
x\$SPWST_FDR=p.adjust(x\$SPWST,method="fdr");
x\$SPWST.up_FDR=p.adjust(x\$SPWST.up,method="fdr");
x\$SPWST.down_FDR=p.adjust(x\$SPWST.down,method="fdr");

write.table(x, "$final_file",quote=FALSE,row.names=FALSE,sep="\t");



out.dir="$input_directory";
if (!file.exists(out.dir)==T) dir.create(out.dir);

q();

_END_OF_R_
    #-------------------------------------------------

    print $fh_R_finisher "$R_command_finisher\n";

    my $cmd = "R --vanilla --slave \< $R_burden_finisher_file";
    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
    );
    unless($return) { 
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }


    return 1;
}


