package Genome::Model::Tools::Validation::CoverageCopynumberPlotWgsVsExome;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatVcf - "Inputs of Varscan and Copynumber, Output of R plots"
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	03-Mar-2011 by W.S.
#	MODIFIED:	24-Apr-2011 by CAM
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Validation::CoverageCopynumberPlotWgsVsExome {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		varscan_file	=> { is => 'Text', doc => "File of varscan validated calls, ex: ", is_optional => 0, is_input => 1 },
		cnvhmm_file	=> { is => 'Text', doc => "File of cnvhmm whole genome predictions", is_optional => 0, is_input => 1 },
		varscan_r_library	=> { is => 'Text', doc => "File of cnvhmm whole genome predictions", is_optional => 0, is_input => 1, default => '/gscmnt/sata423/info/medseq/analysis/CaptureValidationGraphs/VarScanGraphLib.R'},
		sample_id	=> { is => 'Text', doc => "Sample ID to be put on graphs", is_optional => 1, is_input => 1, default => 'unspecified' },
		r_script_output_file     => { is => 'Text', doc => "R script built and run by this module", is_optional => 0, is_input => 1},
		coverage_output_file     => { is => 'Text', doc => "PDF Coverage output file", is_optional => 0, is_input => 1, is_output => 1 },
		copynumber_output_file     => { is => 'Text', doc => "PDF Copynumber output file", is_optional => 0, is_input => 1, is_output => 1 },
		skip_if_output_is_present     => { is => 'Text', doc => "Skip if Output is Present", is_optional => 1, is_input => 1, default => 0},
	],
};

sub sub_command_sort_position { 1 }

sub help_brief {                            # keep this to just a few words <---
    "Inputs of Varscan and Copynumber, Output of R plots"                 
}

sub help_synopsis {
    return <<EOS
Inputs of Varscan and Copynumber, Output of R plots
EXAMPLE:	gmt validation coverage-copynumber-plot --varscan-file input.txt --cnvhmm-file input.cn --copynumber-output-file --coverage-output-file --r-script-output-file --varscan-file --sample-id
EXAMPLE:	gmt validation coverage-copynumber-plot --varscan-file /gscmnt/sata843/info/medseq/wschierd/MMY_Validation/Capture_Validation/MMY1/varscan/SNV_alltiers_manrev.txt --cnvhmm-file /gscmnt/sata843/info/medseq/wschierd/MMY_CopyNumber/MMY1_CNV.seg --copynumber-output-file /gscmnt/sata843/info/medseq/wschierd/MMY_Validation/Capture_Validation/MMY1/varscan/MMY1_SNV_alltiers_Copynumber.pdf --coverage-output-file /gscmnt/sata843/info/medseq/wschierd/MMY_Validation/Capture_Validation/MMY1/varscan/MMY1_SNV_alltiers_Coverage.pdf --r-script-output-file /gscmnt/sata843/info/medseq/wschierd/MMY_Validation/Capture_Validation/MMY1/varscan/R.input --sample-id MMY1

EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	##inputs##
	my $varscan_file = $self->varscan_file;
	my $copynumber_file = $self->cnvhmm_file;
	my $sample_id = $self->sample_id;
	##outputs##
	my $r_script_output_file = $self->r_script_output_file;
	my $coverage_output_file = $self->coverage_output_file;
	my $copynumber_output_file = $self->copynumber_output_file;
	##options##
	my $r_library = $self->varscan_r_library;
        my $skip_if_output_is_present = $self->skip_if_output_is_present;

        my $outfile_varscan_plus_cn;
        if ($varscan_file =~ m/\.txt$/) {
            $outfile_varscan_plus_cn = $varscan_file;
            $outfile_varscan_plus_cn =~ s/\.txt$//;
            $outfile_varscan_plus_cn .= ".copynumber.txt";
        }
        else {
            $outfile_varscan_plus_cn = "$varscan_file".".copynumber";
        }

        if (! -s $outfile_varscan_plus_cn) {
            open(CN_VARSCAN_OUT, ">$outfile_varscan_plus_cn") or die "Can't open output file: $!\n";

            my %copynumber_hash=%{&build_hash($copynumber_file)};

            my $varscan_input = new FileHandle ($varscan_file);
            while (my $line2 = <$varscan_input>) {
                chomp($line2);
                my ($chr, $pos, $ref, $var, $normal_ref, $normal_var, $normal_var_pct, $normal_IUB, $tumor_ref, $tumor_var, $tumor_var_pct, $tumor_IUB, $varscan_call, $germline_pvalue, $somatic_pvalue, @otherstuff) = split(/\t/, $line2);
                my $varscan_cn=&get_cn($chr,$pos,$pos,\%copynumber_hash);
                print CN_VARSCAN_OUT "$line2\t$varscan_cn\n"
            }
            close(CN_VARSCAN_OUT);
        }
        else {
            print STDERR "Found varscan + copy-number file. Continuing...\n";
        }

        # Open Output
        unless (open(R_COMMANDS,">$r_script_output_file")) {
            die "Could not open output file '$r_script_output_file' for writing";
        }
        
        my $readcount_cutoff = 30;
#coverage
#	print R_COMMANDS 'options(echo = FALSE)'."\n";#suppress output to stdout
        print R_COMMANDS 'sink("/dev/null")'."\n";
        print R_COMMANDS "genome=\"$sample_id\";"."\n";
        print R_COMMANDS "source(\"$r_library\");"."\n"; #this contains R functions for loading and graphing VarScan files
        print R_COMMANDS 'library(fpc);'."\n";
        print R_COMMANDS "varscan.load_snp_output(\"$outfile_varscan_plus_cn\",header=F)->xcopy"."\n";
	print R_COMMANDS "varscan.load_snp_output(\"$outfile_varscan_plus_cn\",header=F,min_tumor_depth='$readcount_cutoff',min_normal_depth='$readcount_cutoff')->xcopy2"."\n";
        print R_COMMANDS 'z1=subset(xcopy, xcopy$V13 == "Somatic");'."\n";
	print R_COMMANDS 'z2=subset(xcopy2, xcopy2$V13 == "Somatic");'."\n";
        print R_COMMANDS 'covtum1=(z1$V9+z1$V10);'."\n";
#	print R_COMMANDS 'covtum2=(z2$V9+z2$V10);'."\n";
#	print R_COMMANDS 'maxx=max(c(covtum1,covtum2));'."\n";
        print R_COMMANDS 'covnorm1=(z1$V5+z1$V6);'."\n";
#	print R_COMMANDS 'covnorm2=(z2$V5+z2$V6);'."\n";
#	print R_COMMANDS 'maxx2=max(c(covnorm1,covnorm2));'."\n";
#	print R_COMMANDS 'if (maxx >= 1000) {maxx = 1000};'."\n";
#	print R_COMMANDS 'if (maxx2 >= 1000) {maxx2 = 1000};'."\n";

        #code for making snp inclusion dropoff picture
        print R_COMMANDS 'rc_cutoffs = 1:100; snps_passed_cutoff = NULL;' . "\n";
        print R_COMMANDS 'for (i in rc_cutoffs) { snps_passed_cutoff[i] = dim(z1[(z1$V9+z1$V10) >= i & (z1$V5+z1$V6) >= i,])[1]; }' . "\n";

        #open up image for plotting
        if ($coverage_output_file =~ /.pdf/) {
            print R_COMMANDS "pdf(file=\"$coverage_output_file\",width=10,height=7.5,bg=\"white\");"."\n";
        }
        elsif ($coverage_output_file =~ /.png/) {
            print R_COMMANDS "png(file=\"$coverage_output_file\",width=1200,height=800);"."\n";
        }
        else {
            die "unrecognized coverage output file type...\n";
        }

#	print R_COMMANDS 'par(mfrow=c(2,3));'."\n";
        print R_COMMANDS 'par(mfrow=c(2,3));'."\n";

        #SNP DROPOFF PLOT
        print R_COMMANDS 'plot.default(x=rc_cutoffs,y=snps_passed_cutoff,xlab="Read-count Cut-off",ylab="Number of SNVs", main=paste(genome,"SNVs Passing Read Depth Filter"),cex=0.4);' . "\n";
        print R_COMMANDS "lines(c($readcount_cutoff,$readcount_cutoff),c(0,snps_passed_cutoff[1]), col=\"blue\");"."\n";
        print R_COMMANDS 'legend(x="topright", paste("Filter","Cut-off",sep=" "),col=c("blue"),lty = c(1,1,1), lwd = 1);'."\n"; #top right will rarely have any points

        #TUMOR COVERAGE PLOT
        print R_COMMANDS 'plot.default(x=(z1$V9+z1$V10),y=z1$V11,xlab="WGS Coverage",ylab="Variant Allele Frequency in WGS", main=paste(genome,"WGS SNV Coverage"), type="p",pch=19,cex=0.4,col="#FF000039",xlim=c(0,200),ylim=c(0,100));'."\n";
        print R_COMMANDS 'lines(c(20,20),c(1,100), col="black");'."\n";
        print R_COMMANDS 'lines(c(30,30),c(1,100), col="blue");'."\n";
        print R_COMMANDS 'lines(c(100,100),c(1,100), col="green4");'."\n";
        print R_COMMANDS 'legend(x="topright", title = "Coverage", c("20x", "30x", "100x"),col=c("black","blue","green4"),lty = c(1,1,1), lwd = 1);'."\n"; #top right will rarely have any points

        #NORMAL COVERAGE PLOT
        print R_COMMANDS 'plot.default(x=(z1$V5+z1$V6),y=z1$V7,xlab="Exome Coverage",ylab="Variant Allele Frequency in Exome", main=paste(genome,"Exome SNV Coverage"), type="p",pch=19,cex=0.4,col="#FF000039",xlim=c(0,200),ylim=c(0,100));'."\n";
        print R_COMMANDS 'lines(c(20,20),c(1,100), col="black");'."\n";
        print R_COMMANDS 'lines(c(30,30),c(1,100), col="blue");'."\n";
        print R_COMMANDS 'lines(c(100,100),c(1,100), col="green4");'."\n";
        print R_COMMANDS 'legend(x="topright", title = "Coverage", c("20x", "30x", "100x"),col=c("black","blue","green4"),lty = c(1,1,1), lwd = 1);'."\n"; #top right will rarely have any points

        #TUMOR FREQ WITH CN DENSITY
        #print R_COMMANDS 'plot.default(x=z1$V20,y=z1$V11,xlab="Copy Number",ylab="Variant Allele Frequency in WGS", main=paste(genome,"Copy Number"), type="p",pch=19,cex=0.4,col="#FF000055",xlim=c(0,5),ylim=c(0,100));'."\n";
        #print R_COMMANDS 'plot(density(z1$V20, bw=0.5, from=0,to=5),xlab="Copy Number",ylab="Kernel Density", main=paste(genome,"Copy Number"),lwd=1.5,col="red");'."\n";
        print R_COMMANDS '
cn1minus=subset(z1, z1$V20 >= 0 & z1$V20 <= 1.75);
cn2=subset(z1, z1$V20 >= 1.75 & z1$V20 <= 2.25);
cn3=subset(z1, z1$V20 >= 2.25 & z1$V20 <= 3.5);
cn4plus=subset(z1, z1$V20 >= 3.5);'."\n";

        print R_COMMANDS 'den1 <- 0;'."\n";
	print R_COMMANDS 'den2 <- 0;'."\n";
	print R_COMMANDS 'den3 <- 0;'."\n";
	print R_COMMANDS 'den4 <- 0;'."\n";
        print R_COMMANDS 'den1factor = 0; den2factor = 0; den3factor = 0; den4factor = 0;'."\n";
        print R_COMMANDS 'N = dim(z1)[1];'."\n";
	print R_COMMANDS 'if(dim(cn1minus)[1] < 2) {den1$x = den1$y=1000; den1factor = 0;} else {den1 <- density(cn1minus$V11, from=0,to=100,na.rm=TRUE); den1factor = dim(cn1minus)[1]/N * den1$y;}'."\n";
	print R_COMMANDS 'if(dim(cn2)[1] < 2) {den2$x = den2$y=1000; den2factor = 0;} else {den2 <- density(cn2$V11, from=0,to=100,na.rm=TRUE); den2factor = dim(cn2)[1]/N * den2$y;}'."\n";
	print R_COMMANDS 'if(dim(cn3)[1] < 2) {den3$x = den3$y=1000; den3factor = 0;} else {den3 <- density(cn3$V11, from=0,to=100,na.rm=TRUE); den3factor = dim(cn3)[1]/N * den3$y;}'."\n";
	print R_COMMANDS 'if(dim(cn4plus)[1] < 2) {den4$x = den4$y=1000; den4factor = 0;} else {den4 <- density(cn4plus$V11, from=0,to=100,na.rm=TRUE); den4factor = dim(cn4plus)[1]/N * den4$y;}'."\n";
print R_COMMANDS '
maxden = max(c(den1factor,den2factor,den3factor,den4factor));
finalfactor = 40 / maxden;
#par(mar=c(5,4,8,2) + 0.1);
plot.default(x=cn1minus$V7,y=cn1minus$V11,xlab="Variant Allele Frequency in Exome",ylab="Variant Allele Frequency in WGS", main=paste(genome,"Variant Allele Frequency"), type="p",pch=19,cex=0.4,col="#FF000055",xlim=c(0,100),ylim=c(0,110));
points(x=cn2$V7,y=cn2$V11, type="p",pch=19,cex=0.4,col="#00FF0055");
points(x=cn3$V7,y=cn3$V11, type="p",pch=19,cex=0.4,col="#0000FF55");
points(x=cn4plus$V7,y=cn4plus$V11, type="p",pch=19,cex=0.4,col="#FFA500FF");
lines(c(60,60),c(10,100),lty=2,col="black");
lines((100-(finalfactor * den1factor)),den1$x,col="#FF0000CC",lwd=2);
lines((100-(finalfactor * den2factor)),den2$x,col="#00FF00CC",lwd=2);
lines((100-(finalfactor * den3factor)),den3$x,col="#0000FFCC",lwd=2);
lines((100-(finalfactor * den4factor)),den4$x,col="#FFA500CC",lwd=2);
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(60,100),labels=c(sprintf("%.3f", maxden),0),col="black",tck=0.01);
mtext("Density         ",adj=1, cex=0.7, padj=-0.5);
par(mgp = c(3,1,0));
#par(mar=c(5,4,4,2) + 0.1);
#legend(x="bottomright", title = "Copy Number\nDensity", c("1", "2", "3", "4+"),col=c("#FF0000","#00FF00","#0000FF","#FFA500"),lty = c(1,1,1), lwd = 1);
#legend(x="topleft",horiz=TRUE,xjust=0, c("1", "2", "3", "4+"),col=c("#FF0000","#00FF00","#0000FF","#FFA500"),pch=19,cex=0.9);'."\n";

=cut

text(20,110,labels="0",cex=1);
text(60,110,labels=sprintf("%.3f", maxden),cex=1);
lines(c(20,20),c(107,104),col="black");
lines(c(60,60),c(107,104),col="black");
lines(c(20,60),c(107,107),col="black");
=cut
        #VARIANT FREQUENCY TUMOR VS. NORMAL
        #print R_COMMANDS 'plot.default(x=z1$V7,y=z1$V11,xlab="Exome Variant Allele Frequency",ylab="Variant Allele Frequency in WGS", main=paste(genome,"Allele Frequency"), type="p",pch=19,cex=0.4, col="#FF000039",xlim=c(0,100),ylim=c(0,100));'."\n";


        #VARIANT FREQUENCY TUMOR VS. NORMAL WITH ALL FILTERS
        print R_COMMANDS 'z2_no_CN_events = subset(z2, z2$V20 == 2);' . "\n";
        print R_COMMANDS 'plot.default(x=z2_no_CN_events$V7,y=z2_no_CN_events$V11,xlab="Variant Allele Frequency in Exome",ylab="Variant Allele Frequency in WGS", main=paste(genome,"Allele Frequency With Filters"), type="p",pch=19,cex=0.4, col="#FF000039",xlim=c(0,100),ylim=c(0,100));'."\n";
        print R_COMMANDS "mtext(\"Coverage > $readcount_cutoff, No Copy Number Events\",cex=0.7, padj=-0.5);"."\n";

        #PLOT the density plot
        print R_COMMANDS 'z2_no_CN_events = cbind(z2_no_CN_events,V21=z2_no_CN_events$V11-z2_no_CN_events$V7);'."\n";
        print R_COMMANDS 'plot(density(z2_no_CN_events$V21, from=-20,to=120),col="red",lwd=1.5,xlab="Allele Frequency Difference (WGS - Exome)",ylab="Kernel Density", main=paste(genome,"Allele Frequency Difference"));'."\n";


        print R_COMMANDS "q()\n";

        close R_COMMANDS;

        my $cmd = "R --vanilla --slave \< $r_script_output_file";
        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            output_files => [$coverage_output_file],
            skip_if_output_is_present => $skip_if_output_is_present,
        );
        unless($return) { 
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }
        return $return;
    }

    sub get_cn
    {
        my ($chr,$start,$stop,$hashref)=@_;
        my %info_hash=%{$hashref};
        my $cn;
        foreach my $ch (sort keys %info_hash)
        {
            next unless ($chr eq $ch);
            foreach my $region (sort keys %{$info_hash{$ch}})
            {
                my ($reg_start,$reg_stop)=split/\_/,$region;
                if ($reg_start<=$start && $reg_stop>=$stop)
                {
                    $cn=$info_hash{$ch}{$region};
                    last;

                }
            }
        }

        $cn=2 unless ($cn);
        return $cn;
    }


    sub build_hash
    {
        my ($file)=@_;
        my %info_hash;
        my $fh=new FileHandle($file);
        while(my $line = <$fh>)
        {
            chomp($line);
            unless ($line =~ /^\w+\t\d+\t\d+\t/) { next;}
            my ($chr,$start,$end,$size,$nmarkers,$tcn,$adjusted_tcn,$ncn,$adjusted_ncn,$score,$status)=split(/\t/, $line);
            if ($score >= 100) {
                my $pos=$start."_".$end;
                $info_hash{$chr}{$pos}=$adjusted_tcn;
            }
        }
        return \%info_hash;
    }
