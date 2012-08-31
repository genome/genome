package Genome::Model::Tools::Validation::CoverageCopynumberPlot;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatVcf - "Inputs of Varscan and Copynumber, Output of R plots"
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	03-Mar-2011 by W.S.
#	MODIFIED:	03-Mar-2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Validation::CoverageCopynumberPlot {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		varscan_file	=> { is => 'Text', doc => "File of varscan validated calls, ex: ", is_optional => 0, is_input => 1 },
		cnvhmm_file	=> { is => 'Text', doc => "File of cnvhmm whole genome predictions", is_optional => 0, is_input => 1 },
		varscan_r_library	=> { is => 'Text', doc => "File of cnvhmm whole genome predictions", is_optional => 0, is_input => 1, default => '/gscmnt/sata423/info/medseq/analysis/CaptureValidationGraphs/VarScanGraphLib.R'},
		sample_id	=> { is => 'Text', doc => "Sample ID to be put on graphs", is_optional => 1, is_input => 1, default => 'unspecified' },
		analysis_type      => { is => 'Text', doc => "Either \'wgs\' for somatic pipeline output or \'capture\' for validation pipeline output", is_optional => 1, is_input => 1, default => 'capture'},
#		readcount_cutoff      => { is => 'Text', doc => "Choose a depth cutoff, this should be 100 for validation and 30 for wgs", is_optional => 1, is_input => 1, default => '100'},
		chr_highlight      => { is => 'Text', doc => "Choose a Chromosome to Highlight with Purple Circles on Plot", is_optional => 1, is_input => 1, default => 'X'},
		positions_highlight      => { is => 'Text', doc => "A tab-delim file list of positions chr\\tposition to highlight on plots", is_optional => 1, is_input => 1},
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
	This tool is used to assess the sources of variation of variant allele freqency within the tumor sample. This is done through an analysis involving both coverage and copy number, where underlying allele frequency is estimated to be reflected as the percentage of reads aligned at a particular location that report a specific base. This analysis has a significant impact on determining whether those sources of variation come from purity, clonality, or various other unexplained factors. The assumption made here is that most CN neutral events with proper coverage should be heterozygous, expressed by reporting half of the reads as the reference allele and half of the reads as the variant allele. In a few cases, the variant mutation will be homozygous, causing all aligned reads to report the allele as variant. Any deviation from these expectations come from factors beyond the events described in the tool, but that can be further investigated. The method used to assess these variations is through assessment of peaks on a kernal density plot. Kernel density estimation (KDE) is used to smooth the data in a non-parametric way to estimate the probability density function of the allele frequencies. In this context, since true allelic frequencies are assumed to be under an error model that causes a normal distribution of allelic freqencies around the true frequency, peaks in KDE plots are used to represent allelic freqencies that are highly reported in the somatic events, and therefore are believed to be the true underlying variant allele frequency.
	In order to have sufficient data (large n), the data at any particular location must achieve adequate coverage. This is because the error in allelic frequency estimation is under strong forces of bias at low coverage. In seqence data it is very likely that due to factors such as GC content and repeptitive regions, that certain chromosomal regions fail to achieve more than a low-depth of coverage. In those cases, it is hard to assess the variables impacting variant allele frequency, as the sequencing error rate becomes greater than any of these subtle effects for which we are trying to measure. In testing, it has been shown that in high-depth data (that above 100x coverage) the variant allele freqency differences are no longer significantly biased by sequencing variability. This is why we use 100x as a cutoff for analysis, with sequence data required to meet or exceed this threshold in a large majority of targeted positions in order for a sample to be considered for this analysis (SEE FIGURE XYZ, THE SIGMOIDAL PLOT OF LOSS OF SNVS AT 100X CUTOFF).
	Once proper coverage is attained, one can assess the causes of variations from expected variant allele frequency. The first possible cause is copy number. Each variant will have frequency distributions according to the number of mutant alleles divided by the total number of alleles. It is for this reason that regions that have a CN alteration will fail to have 50% variant allele frequency, but instead have a distribution determined by the resulting CN. For example, regions with a CN of 3 would be expected to have possible allele frequencies of 33, 66, and 100%.
	Another cause of variations in variant allele frequency is cellular heterogeneity within the tumor sample. This happens in one of two ways: 1) If the mutations are homozygous but has a contamination factor of normal cells within the tumor cell sample, or 2) when there are heterozygous mutations but the tumor cell population has sub-cellularity (clonality) that do not have all of the somatic mutations that are being reported. In the first case, purity can be a significant cause of sampling bias in those tumor types that cannot accurately determine normal versus tumor cells. Most solid tumors will experience effects from this problem, as it can be difficult to be precise in cellular selection. Liquid tumors rarely will experience this bias in the tumor cell population, as they are generally much easier to sort-select. In the second case, tumor clonality has long been a question unanswered by sequencing technologies. In order to accurately assess clonality, there must be deep coverage (100x or more) at a locus, and that locus must have an accurately determined CN status. Once that is achieved, this tool plots the kernal density of the events to look for peaks that are significantly altered from the expected variant allele frequency for any given CN status.
	Possible sources of bias to this examination are any large-scale CN neutral event such as CN-neutral LOH. In this case, even a pure, homogeneous selection of cells can lead to variant allele frequencies that differ from the expected. The hypothesis for this sort of event would be that the somatic event has arisen after an LOH event. This can be assessed through an analysis of allelic imbalance (when the tumor sample allele frequencies are near 0 or 100% while the normal cell sampling reports an allele frequency at or near 50%). In this case, if the expected allele freqencies are far from the reported allele freqencies, it can be attributed to this bias if there is a significant source of allelic imbalance, but that somatic mutation around these events are still in less than 50% of the reads.
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
	my $readcount_cutoff;
#       my $readcount_cutoff = $self->readcount_cutoff;
	my $chr_highlight = $self->chr_highlight;
	my $positions_highlight = $self->positions_highlight;
	##outputs##
	my $r_script_output_file = $self->r_script_output_file;
	my $coverage_output_file = $self->coverage_output_file;
	my $copynumber_output_file = $self->copynumber_output_file;
	##options##
	my $r_library = $self->varscan_r_library;
	my $skip_if_output_is_present = $self->skip_if_output_is_present;
	my $analysis_type = $self->analysis_type;
	if ($analysis_type eq 'wgs') {
		$readcount_cutoff = 20;
	}
	elsif ($analysis_type eq 'capture') {
		$readcount_cutoff = 100;
	}
	else {
		die "analysis type: $analysis_type not supported, choose either wgs or capture";
	}
	my $position_added = 0;
	my %position_highlight_hash;
	if ($positions_highlight && -s $positions_highlight) {
		my $positions_input = new FileHandle ($positions_highlight);
		while (my $line2 = <$positions_input>) {
			$position_added++;
			chomp($line2);
			my ($chr, $pos) = split(/\t/, $line2);
			my $matcher = "$chr\t$pos";
			$position_highlight_hash{$matcher}++;
		}
	}

	## Build temp file for positions where readcounts are needed ##
	my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
	unless($tfh) {
		$self->error_message("Unable to create temporary file $!");
		die;
	}
	$temp_path =~ s/\:/\\\:/g;

	## Build temp file for extra positions to highlight ##
	my ($tfh2,$temp_path2) = Genome::Sys->create_temp_file;
	unless($tfh2) {
		$self->error_message("Unable to create temporary file $!");
		die;
	}
	$temp_path2 =~ s/\:/\\\:/g;

	my %copynumber_hash_tumor=%{&build_hash($copynumber_file,'tumor')};
	my %copynumber_hash_normal=%{&build_hash($copynumber_file,'normal')};
	my $varscan_input = new FileHandle ($varscan_file);
	while (my $line2 = <$varscan_input>) {
		chomp($line2);
		my ($chr, $pos, $ref, $var, $normal_ref, $normal_var, $normal_var_pct, $normal_IUB, $tumor_ref, $tumor_var, $tumor_var_pct, $tumor_IUB, $varscan_call, $germline_pvalue, $somatic_pvalue, @otherstuff) = split(/\t/, $line2);
		my $varscan_cn_tumor=&get_cn($chr,$pos,$pos,\%copynumber_hash_tumor);
		my $varscan_cn_normal=&get_cn($chr,$pos,$pos,\%copynumber_hash_normal);
		print $tfh "$line2\t$varscan_cn_tumor\t$varscan_cn_normal\n";
		if ($positions_highlight && -s $positions_highlight) {
			my $matcher = "$chr\t$pos";
			if (defined $position_highlight_hash{$matcher}) {
				$position_added--;
				my $depth = $tumor_ref + $tumor_var;
				my $varallelefreq = $tumor_var_pct;
				$varallelefreq =~ s/%//;
				print $tfh2 "$varscan_cn_tumor\t$varallelefreq\t$depth\n";
			}
		}
	}
	$tfh->close;
	$tfh2->close;

	if ($positions_highlight) {
		unless($position_added == 0) {
			warn "There are positions in positions_highlight file that aren't present in varscan file...check files for accuracy for missing positions";
		}
	}

	# Open Output
	unless (open(R_COMMANDS,">$r_script_output_file")) {
	    die "Could not open output file '$r_script_output_file' for writing";
	  }

#coverage
	#set defined cutoffs for graphs
	my $maxx = my $absmaxx = 0;
	if ($analysis_type eq 'wgs') {
		$maxx = 200;
		$absmaxx = 500;
	}
	elsif ($analysis_type eq 'capture') {
		$maxx = 1000;
		$absmaxx = 5000;
	}

        my $R_command = <<"_END_OF_R_";
#options(echo = FALSE);#suppress output to stdout
#sink("/dev/null");
genome=\"$sample_id\";
source(\"$r_library\"); #this contains R functions for loading and graphing VarScan files
library(fpc);
library(scatterplot3d);
varscan.load_snp_output(\"$temp_path\",header=F)->xcopy;
varscan.load_snp_output(\"$temp_path\",header=F,min_tumor_depth=$readcount_cutoff,min_normal_depth=$readcount_cutoff)->xcopy100;

additional_plot_points = 0;
_END_OF_R_
        print R_COMMANDS "$R_command\n";
	if ($positions_highlight && -s $positions_highlight && 1 && 1) { #these &&1 mean nothing, they just make my text editor color things correctly (it hates -s without being s///)
$R_command = <<"_END_OF_R_";
additional_plot_points <- read.table(\"$temp_path2\", header = FALSE, sep = "\t");
additional_plot_points_cn1=subset(additional_plot_points, additional_plot_points\$V1 >= 0 & additional_plot_points\$V1 <= 1.75);
additional_plot_points_cn2=subset(additional_plot_points, additional_plot_points\$V1 >= 1.75 & additional_plot_points\$V1 <= 2.25);
additional_plot_points_cn3=subset(additional_plot_points, additional_plot_points\$V1 >= 2.25 & additional_plot_points\$V1 <= 3.5);
additional_plot_points_cn4=subset(additional_plot_points, additional_plot_points\$V1 >= 3.5);
_END_OF_R_
        print R_COMMANDS "$R_command\n";
	}
        $R_command = <<"_END_OF_R_";

z1=subset(xcopy, xcopy\$V13 == "Somatic");
z2=subset(xcopy100, xcopy100\$V13 == "Somatic");
xchr=subset(z1,z1\$V1 == "$chr_highlight");
xchr100=subset(z2,z2\$V1 == "$chr_highlight");
covtum1=(z1\$V9+z1\$V10);
covtum2=(z2\$V9+z2\$V10);
absmaxx=maxx=max(c(covtum1,covtum2));
covnorm1=(z1\$V5+z1\$V6);
covnorm2=(z2\$V5+z2\$V6);
absmaxx2=maxx2=max(c(covnorm1,covnorm2));

#if (maxx >= 1200) {maxx = 1200};
#if (maxx2 >= 1200) {maxx2 = 1200};
#if (maxx <= 800) {maxx = 800};
#if (maxx2 <= 800) {maxx2 = 800};
#if (absmaxx <= 5000) {absmaxx = 5000};
#if (absmaxx2 <= 5000) {absmaxx2 = 5000};

maxx = $maxx;
maxx2 = $maxx;
absmaxx = $absmaxx;
absmaxx2 = $absmaxx;

cn1minus=subset(z1, z1\$V20 >= 0 & z1\$V20 <= 1.75);
cn2=subset(z1, z1\$V20 >= 1.75 & z1\$V20 <= 2.25);
cn3=subset(z1, z1\$V20 >= 2.25 & z1\$V20 <= 3.5);
cn4plus=subset(z1, z1\$V20 >= 3.5);
cn1minus100x=subset(z2, z2\$V20 >= 0 & z2\$V20 <= 1.75);
cn2100x=subset(z2, z2\$V20 >= 1.75 & z2\$V20 <= 2.25);
cn3100x=subset(z2, z2\$V20 >= 2.25 & z2\$V20 <= 3.5);
cn4plus100x=subset(z2, z2\$V20 >= 3.5);

cn1xchr=subset(xchr, xchr\$V20 >= 0 & xchr\$V20 <= 1.75);
cn2xchr=subset(xchr, xchr\$V20 >= 1.75 & xchr\$V20 <= 2.25);
cn3xchr=subset(xchr, xchr\$V20 >= 2.25 & xchr\$V20 <= 3.5);
cn4xchr=subset(xchr, xchr\$V20 >= 3.5);
cn1xchr100=subset(xchr100, xchr100\$V20 >= 0 & xchr100\$V20 <= 1.75);
cn2xchr100=subset(xchr100, xchr100\$V20 >= 1.75 & xchr100\$V20 <= 2.25);
cn3xchr100=subset(xchr100, xchr100\$V20 >= 2.25 & xchr100\$V20 <= 3.5);
cn4xchr100=subset(xchr100, xchr100\$V20 >= 3.5);

cov20x=subset(z1, (z1\$V9+z1\$V10) <= 20);
cov50x=subset(z1, (z1\$V9+z1\$V10) >= 20 & (z1\$V9+z1\$V10) <= 50);
cov100x=subset(z1, (z1\$V9+z1\$V10) >= 50 & (z1\$V9+z1\$V10) <= 100);
cov100xplus=subset(z1, (z1\$V9+z1\$V10) >= 100);

den1 <- 0;
den2 <- 0;
den3 <- 0;
den4 <- 0;
den1100x <-  0;
den2100x <-  0;
den3100x <-  0;
den4100x <-  0;

den1factor = 0; den2factor = 0; den3factor = 0; den4factor = 0;
den1factor100 = 0; den2factor100 = 0; den3factor100 = 0; den4factor100 = 0;

N = dim(z1)[1];
N100 = dim(z2)[1];

if(dim(cn1minus)[1] < 2) {den1\$x = den1\$y=1000;} else {den1 <- density(cn1minus\$V11, from=0,to=100,na.rm=TRUE); den1factor = dim(cn1minus)[1]/N * den1\$y;};
if(dim(cn2)[1] < 2) {den2\$x = den2\$y=1000;} else {den2 <- density(cn2\$V11, from=0,to=100,na.rm=TRUE); den2factor = dim(cn2)[1]/N * den2\$y;};
if(dim(cn3)[1] < 2) {den3\$x = den3\$y=1000;} else {den3 <- density(cn3\$V11, from=0,to=100,na.rm=TRUE); den3factor = dim(cn3)[1]/N * den3\$y;};
if(dim(cn4plus)[1] < 2) {den4\$x = den4\$y=1000;} else {den4 <- density(cn4plus\$V11, from=0,to=100,na.rm=TRUE); den4factor = dim(cn4plus)[1]/N * den4\$y;};
if(dim(cn1minus100x)[1] < 2) {den1100x\$x = den1100x\$y=1000} else {den1100x <- density(cn1minus100x\$V11, from=0,to=100,na.rm=TRUE); den1factor100 = dim(cn1minus100x)[1]/N100 * den1100x\$y;};
if(dim(cn2100x)[1] < 2) {den2100x\$x = den2100x\$y=1000} else {den2100x <- density(cn2100x\$V11, from=0,to=100,na.rm=TRUE);den2factor100 = dim(cn2100x)[1]/N100 * den2100x\$y;};
if(dim(cn3100x)[1] < 2) {den3100x\$x = den3100x\$y=1000} else {den3100x <- density(cn3100x\$V11, from=0,to=100,na.rm=TRUE);den3factor100 = dim(cn3100x)[1]/N100 * den3100x\$y;};
if(dim(cn4plus100x)[1] < 2) {den4100x\$x = den4100x\$y=1000} else {den4100x <- density(cn4plus100x\$V11, from=0,to=100,na.rm=TRUE);den4factor100 = dim(cn4plus100x)[1]/N100 * den4100x\$y;};

dennormcov <- density((z1\$V5+z1\$V6), bw=4, from=0,to=maxx,na.rm=TRUE);
dentumcov <- density((z1\$V9+z1\$V10), bw=4, from=0,to=maxx,na.rm=TRUE);
dennormcov100x <- density((z2\$V5+z2\$V6), bw=4, from=0,to=maxx,na.rm=TRUE);
dentumcov100x <- density((z2\$V9+z2\$V10), bw=4, from=0,to=maxx,na.rm=TRUE);

#find inflection points (peaks)
#den2diff = diff(den2\$y);

peaks<-function(series,span=3)
{
z <- embed(series, span);
s <- span%/%2;
v<- max.col(z) == 1 + s;
result <- c(rep(FALSE,s),v);
result <- result[1:(length(result)-s)];
result;
} 

#labels to use for density plot values
if(dim(cn1minus)[1] < 2) {cn1peaks = cn1peakpos = cn1peakheight = 0;} else {cn1peaks = peaks(den1factor); cn1peaks = append(cn1peaks,c("FALSE","FALSE"),after=length(cn1peaks)); cn1peakpos = subset(den1\$x,cn1peaks==TRUE & den1\$y > 0.001); cn1peakheight = subset(den1factor,cn1peaks==TRUE & den1\$y > 0.001);}
if(dim(cn2)[1] < 2) {cn2peaks = cn2peakpos = cn2peakheight = 0;} else {cn2peaks = peaks(den2factor); cn2peaks = append(cn2peaks,c("FALSE","FALSE"),after=length(cn2peaks)); cn2peakpos = subset(den2\$x,cn2peaks==TRUE & den2\$y > 0.001); cn2peakheight = subset(den2factor,cn2peaks==TRUE & den2\$y > 0.001);}
if(dim(cn3)[1] < 2) {cn3peaks = cn3peakpos = cn3peakheight = 0;} else {cn3peaks = peaks(den3factor); cn3peaks = append(cn3peaks,c("FALSE","FALSE"),after=length(cn3peaks)); cn3peakpos = subset(den3\$x,cn3peaks==TRUE & den3\$y > 0.001); cn3peakheight = subset(den3factor,cn3peaks==TRUE & den3\$y > 0.001);}
if(dim(cn4plus)[1] < 2) {cn4peaks = cn4peakpos = cn4peakheight = 0;} else {cn4peaks = peaks(den4factor); cn4peaks = append(cn4peaks,c("FALSE","FALSE"),after=length(cn4peaks)); cn4peakpos = subset(den4\$x,cn4peaks==TRUE & den4\$y > 0.001); cn4peakheight = subset(den4factor,cn4peaks==TRUE & den4\$y > 0.001);}
if(dim(cn1minus100x)[1] < 2) {cn1peaks100 = cn1peakpos100 = cn1peakheight100 = 0;} else {cn1peaks100 = peaks(den1factor100); cn1peaks100 = append(cn1peaks100,c("FALSE","FALSE"),after=length(cn1peaks100)); cn1peakpos100 = subset(den1100x\$x,cn1peaks100==TRUE & den1100x\$y > 0.001); cn1peakheight100 = subset(den1factor100,cn1peaks100==TRUE & den1100x\$y > 0.001);}
if(dim(cn2100x)[1] < 2) {cn2peaks100 = cn2peakpos100 = cn2peakheight100 = 0;} else {cn2peaks100 = peaks(den2factor100); cn2peaks100 = append(cn2peaks100,c("FALSE","FALSE"),after=length(cn2peaks100)); cn2peakpos100 = subset(den2100x\$x,cn2peaks100==TRUE & den2100x\$y > 0.001); cn2peakheight100 = subset(den2factor100,cn2peaks100==TRUE & den2100x\$y > 0.001);}
if(dim(cn3100x)[1] < 2) {cn3peaks100 = cn3peakpos100 = cn3peakheight100 = 0;} else {cn3peaks100 = peaks(den3factor100); cn3peaks100 = append(cn3peaks100,c("FALSE","FALSE"),after=length(cn3peaks100)); cn3peakpos100 = subset(den3100x\$x,cn3peaks100==TRUE & den3100x\$y > 0.001); cn3peakheight100 = subset(den3factor100,cn3peaks100==TRUE & den3100x\$y > 0.001);}
if(dim(cn4plus100x)[1] < 2) {cn4peaks100 = cn4peakpos100 = cn4peakheight100 = 0;} else {cn4peaks100 = peaks(den4factor100); cn4peaks100 = append(cn4peaks100,c("FALSE","FALSE"),after=length(cn4peaks100)); cn4peakpos100 = subset(den4100x\$x,cn4peaks100==TRUE & den4100x\$y > 0.001); cn4peakheight100 = subset(den4factor100,cn4peaks100==TRUE & den4100x\$y > 0.001);}

_END_OF_R_

        print R_COMMANDS "$R_command\n";


        #open up image for plotting
        if ($coverage_output_file =~ /.pdf/) {
            print R_COMMANDS "pdf(file=\"$coverage_output_file\",width=10,height=7.5,bg=\"white\");"."\n";
        }
        elsif ($coverage_output_file =~ /.png/) {
            print R_COMMANDS "png(file=\"$coverage_output_file\",width=1200,height=800);"."\n";
        }
        else {
            die "unrecognized coverage output file type...please append .pdf or .png to the end of your coverage output file\n";
        }

        $R_command = <<"_END_OF_R_";

par(mfrow=c(2,3));

        #NORMAL COVERAGE PLOT
plot.default(x=(cn1minus\$V5+cn1minus\$V6),y=(cn1minus\$V7),xlab="Normal Coverage",ylab="Normal Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#1C366044",xlim=c(0,maxx2),ylim=c(0,100));
points(x=(cn2\$V5+cn2\$V6),y=(cn2\$V7), type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=(cn3\$V5+cn3\$V6),y=(cn3\$V7), type="p",pch=19,cex=0.4,col="#F4981955");
points(x=(cn4plus\$V5+cn4plus\$V6),y=(cn4plus\$V7), type="p",pch=19,cex=0.4,col="#E5242044");
points(dennormcov\$x,((dennormcov\$y * 1000)+20),col="#0000000F", type="p",pch=19,cex=0.4);
#label x
points(x=(xchr\$V5+xchr\$V6),y=(xchr\$V7), type="p",pch=1,cex=0.8,col="#A020F099");
lines(c(20,20),c(1,100), col="black");
lines(c(30,30),c(1,100), col="blue");
lines(c(100,100),c(1,100), col="green4");
legend(x="topright", title = "Copy Number", xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,1),cex=0.6);
legend(x="right", title = "Coverage", c("20x", "30x", "100x", "N"),col=c("black","blue","green4","#00000055"),lty = c(1,1,1,1), lwd = c(1,1,1,2),cex=0.6);

        #TUMOR CN PLOT
maxden = max(c(den1factor,den2factor,den3factor,den4factor));
finalfactor = 40 / maxden;
plot.default(x=cn1minus\$V7,y=cn1minus\$V11,xlab="Variant Allele Frequency in Normal",ylab="Variant Allele Frequency in Tumor", main=paste(genome,"Variant Allele Frequency"), type="p",pch=19,cex=0.4,col="#1C366044",xlim=c(0,100),ylim=c(0,110));
points(x=cn2\$V7,y=cn2\$V11, type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=cn3\$V7,y=cn3\$V11, type="p",pch=19,cex=0.4,col="#F4981955");
points(x=cn4plus\$V7,y=cn4plus\$V11, type="p",pch=19,cex=0.4,col="#E5242055");
#label x
points(x=xchr\$V7,y=xchr\$V11, type="p",pch=1,cex=0.8,col="#A020F099");
lines(c(60,60),c(10,100),lty=2,col="black");
lines((100-(finalfactor * den1factor)),den1\$x,col="#1C3660AA",lwd=2);
lines((100-(finalfactor * den2factor)),den2\$x,col="#67B32EAA",lwd=2);
lines((100-(finalfactor * den3factor)),den3\$x,col="#F49819AA",lwd=2);
lines((100-(finalfactor * den4factor)),den4\$x,col="#E52420AA",lwd=2);
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(60,100),labels=c(sprintf("%.3f", maxden),0),col="black",tck=0.01);
mtext("CN Density         ",adj=1, cex=0.7, padj=-0.5);
par(mgp = c(3,1,0));
#par(mar=c(5,4,4,2) + 0.1);
legend(x="topleft",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,1),cex=0.6);

        #TUMOR COVERAGE PLOT
#plotting coverage and copynumber on log scale with density plot incorporated with x chr labels
finalfactor = 25 / maxden;
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
points(x=(cn2\$V9+cn2\$V10),y=(cn2\$V11),type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=(cn1minus\$V9+cn1minus\$V10),y=(cn1minus\$V11),type="p",pch=19,cex=0.4,col="#1C366044");
points(x=(cn3\$V9+cn3\$V10),y=(cn3\$V11),type="p",pch=19,cex=0.4,col="#F4981955");
points(x=(cn4plus\$V9+cn4plus\$V10),y=(cn4plus\$V11),type="p",pch=19,cex=0.4,col="#E5242044");
#label x
points(x=(cn2xchr\$V9+cn2xchr\$V10),y=(cn2xchr\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
points(x=(cn1xchr\$V9+cn1xchr\$V10),y=(cn1xchr\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
points(x=(cn3xchr\$V9+cn3xchr\$V10),y=(cn3xchr\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
points(x=(cn4xchr\$V9+cn4xchr\$V10),y=(cn4xchr\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
par(new=TRUE);
plot.default(x=c(1:10),y=c(1:10),ylim=c(0,110),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000");
lines(c(25,25),c(10,100),lty=2,col="black");
lines(((finalfactor * den2factor)),den2\$x,col="#67B32EAA",lwd=2);
lines(((finalfactor * den1factor)),den1\$x,col="#1C3660AA",lwd=2);
lines(((finalfactor * den3factor)),den3\$x,col="#F49819AA",lwd=2);
lines(((finalfactor * den4factor)),den4\$x,col="#E52420AA",lwd=2);
text(x=(finalfactor * cn1peakheight)+2,y=cn1peakpos,labels=signif(cn1peakpos,3),cex=0.7,srt=-90,col="#1C3660AA");
text(x=(finalfactor * cn3peakheight)+2,y=cn3peakpos,labels=signif(cn3peakpos,3),cex=0.7,srt=-90,col="#F49819AA");
text(x=(finalfactor * cn4peakheight)+2,y=cn4peakpos,labels=signif(cn4peakpos,3),cex=0.7,srt=-90,col="#E52420AA");
text(x=(finalfactor * cn2peakheight)+2,y=cn2peakpos,labels=signif(cn2peakpos,3),cex=0.7,srt=-90,col="#67B32EAA");
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(0,25),labels=c(0,sprintf("%.3f", maxden)),col="black",tck=0.01);
mtext("         CN Density",adj=0, cex=0.7, padj=-0.5);
par(mgp = c(3,1,0));
legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);

        #SNP DROPOFF PLOT
        #code for making snp inclusion dropoff picture
rc_cutoffs = 1:800; snps_passed_cutoff = NULL;
for (i in rc_cutoffs) { snps_passed_cutoff[i] = dim(z1[(z1\$V9+z1\$V10) >= i & (z1\$V5+z1\$V6) >= i,])[1]; }
plot.default(x=rc_cutoffs,y=snps_passed_cutoff,xlab="Read-count Cut-off",ylab="Number of SNVs", main=paste(genome,"SNVs Passing Read Depth Filter"),cex=0.4);
lines(c($readcount_cutoff,$readcount_cutoff),c(0,snps_passed_cutoff[1]), col=\"blue\");
legend(x="topright", paste("Filter","Cut-off",sep=" "),col=c("blue"),lty = c(1,1,1), lwd = 1);

        #BEGIN COVERAGE > 100 PLOTS
        #TUMOR CN PLOT 100x
maxden100 = max(c(den1factor100,den2factor100,den3factor100,den4factor100));
finalfactor = 40 / maxden100;
plot.default(x=cn1minus100x\$V7,y=cn1minus100x\$V11,xlab="Normal Variant Allele Frequency",ylab="Tumor Variant Allele Frequency", main=paste(genome," Allele Frequency"),type="p",pch=19,cex=0.4,col="#1C366044",xlim=c(0,100),ylim=c(0,110));
points(x=cn2100x\$V7,y=cn2100x\$V11, type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=cn3100x\$V7,y=cn3100x\$V11, type="p",pch=19,cex=0.4,col="#F4981955");
points(x=cn4plus100x\$V7,y=cn4plus100x\$V11, type="p",pch=19,cex=0.4,col="#E5242044");
#label x
points(x=xchr100\$V7,y=xchr100\$V11, type="p",pch=1,cex=0.8,col="#A020F099");
lines(c(60,60),c(10,100),lty=2,col="black");
lines((100-(finalfactor * den1factor100)),den1100x\$x,col="#1C3660AA",lwd=2);
lines((100-(finalfactor * den2factor100)),den2100x\$x,col="#67B32EAA",lwd=2);
lines((100-(finalfactor * den3factor100)),den3100x\$x,col="#F49819AA",lwd=2);
lines((100-(finalfactor * den4factor100)),den4100x\$x,col="#E52420AA",lwd=2);
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(60,100),labels=c(sprintf("%.3f", maxden100),0),col="black",tck=0.01);
#mtext("CN Density         ",adj=1, cex=0.7, padj=-0.5);
par(mgp = c(3,1,0));
#par(mar=c(5,4,4,2) + 0.1);
legend(x="topleft",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,1),cex=0.6);
mtext("Normal and Tumor Coverage > $readcount_cutoff",cex=0.7, padj=-0.5);

        #TUMOR COVERAGE PLOT 100x
#plotting coverage and copynumber on log scale with density plot incorporated with x chr labels
finalfactor = 25 / maxden100;
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
points(x=(cn2100x\$V9+cn2100x\$V10),y=(cn2100x\$V11),type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=(cn1minus100x\$V9+cn1minus100x\$V10),y=(cn1minus100x\$V11),type="p",pch=19,cex=0.4,col="#1C366044");
points(x=(cn3100x\$V9+cn3100x\$V10),y=(cn3100x\$V11),type="p",pch=19,cex=0.4,col="#F4981999");
points(x=(cn4plus100x\$V9+cn4plus100x\$V10),y=(cn4plus100x\$V11),type="p",pch=19,cex=0.4,col="#E5242044");
#label x
points(x=(cn2xchr100\$V9+cn2xchr100\$V10),y=(cn2xchr100\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
points(x=(cn1xchr100\$V9+cn1xchr100\$V10),y=(cn1xchr100\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
points(x=(cn3xchr100\$V9+cn3xchr100\$V10),y=(cn3xchr100\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
points(x=(cn4xchr100\$V9+cn4xchr100\$V10),y=(cn4xchr100\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
par(new=TRUE);
plot.default(x=c(1:10),y=c(1:10),ylim=c(0,110),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000");
lines(c(25,25),c(10,100),lty=2,col="black");
lines(((finalfactor * den2factor100)),den2100x\$x,col="#67B32EAA",lwd=2);
lines(((finalfactor * den1factor100)),den1100x\$x,col="#1C3660AA",lwd=2);
lines(((finalfactor * den3factor100)),den3100x\$x,col="#F49819AA",lwd=2);
lines(((finalfactor * den4factor100)),den4100x\$x,col="#E52420AA",lwd=2);
text(x=(finalfactor * cn1peakheight100)+2,y=cn1peakpos100,labels=signif(cn1peakpos100,3),cex=0.7,srt=-90,col="#1C3660AA");
text(x=(finalfactor * cn3peakheight100)+2,y=cn3peakpos100,labels=signif(cn3peakpos100,3),cex=0.7,srt=-90,col="#F49819AA");
text(x=(finalfactor * cn4peakheight100)+2,y=cn4peakpos100,labels=signif(cn4peakpos100,3),cex=0.7,srt=-90,col="#E52420AA");
text(x=(finalfactor * cn2peakheight100)+2,y=cn2peakpos100,labels=signif(cn2peakpos100,3),cex=0.7,srt=-90,col="#67B32EAA");
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(0,25),labels=c(0,sprintf("%.3f", maxden100)),col="black",tck=0.01);
par(mgp = c(3,1,0));
legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);
mtext("Normal and Tumor Coverage > $readcount_cutoff",cex=0.7, padj=-0.5);

        #TUMOR COVERAGE PLOT, OLD STYLE VARIATIONS
plot.default(x=(cn1minus\$V9+cn1minus\$V10),y=(cn1minus\$V11),xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#1C366044",xlim=c(0,maxx),ylim=c(0,110));
points(x=(cn2\$V9+cn2\$V10),y=(cn2\$V11),type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=(cn3\$V9+cn3\$V10),y=(cn3\$V11),type="p",pch=19,cex=0.4,col="#F4981955");
points(x=(cn4plus\$V9+cn4plus\$V10),y=(cn4plus\$V11),type="p",pch=19,cex=0.4,col="#E5242044");
points(dentumcov\$x,((dentumcov\$y * 1000)),col="#0000000F", type="p",pch=19,cex=0.4);
#label x
points(x=(xchr\$V9+xchr\$V10),y=(xchr\$V11),type="p",pch=1,cex=0.8,col="#A020F099");
lines(c(20,20),c(1,100), col="black");
lines(c(30,30),c(1,100), col="blue");
lines(c(100,100),c(1,100), col="green4");
legend(x="topleft",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,1),cex=0.6);
legend(x="topright", title = "Coverage", c("20x", "30x", "100x", "N"),col=c("black","blue","green4","#00000055"),lty = c(1,1,1,1), lwd = c(1,1,1,2),cex=0.6);

#100x coverage
plot.default(x=(cn1minus100x\$V9+cn1minus100x\$V10),y=(cn1minus100x\$V11),xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#1C366044",xlim=c(0,maxx),ylim=c(0,110));
points(x=(cn2100x\$V9+cn2100x\$V10),y=(cn2100x\$V11), type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=(cn3100x\$V9+cn3100x\$V10),y=(cn3100x\$V11), type="p",pch=19,cex=0.4,col="#F4981955");
points(x=(cn4plus100x\$V9+cn4plus100x\$V10),y=(cn4plus100x\$V11), type="p",pch=19,cex=0.4,col="#E5242044");
points(dentumcov100x\$x,((dentumcov100x\$y * 1000)),col="#0000000F", type="p",pch=19,cex=0.4);
#label x
points(x=(xchr100\$V9+xchr100\$V10),y=(xchr100\$V11),type="p",pch=1,cex=0.8,col="#A020F099");
mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);
lines(c(20,20),c(1,100), col="black");
lines(c(30,30),c(1,100), col="blue");
lines(c(100,100),c(1,100), col="green4");
legend(x="topleft",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,1),cex=0.6);
legend(x="topright", title = "Coverage", c("20x", "30x", "100x", "N"),col=c("black","blue","green4","#00000055"),lty = c(1,1,1,1), lwd = c(1,1,1,2),cex=0.6);


##########################################################

#transparency changes
#plotting coverage and copynumber on log scale with density plot incorporated with x chr labels
finalfactor = 25 / maxden;
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
points(x=(cn2\$V9+cn2\$V10),y=(cn2\$V11),type="p",pch=19,cex=0.4,col="#67B32E");
points(x=(cn1minus\$V9+cn1minus\$V10),y=(cn1minus\$V11),type="p",pch=19,cex=0.4,col="#1C3660");
points(x=(cn3\$V9+cn3\$V10),y=(cn3\$V11),type="p",pch=19,cex=0.4,col="#F49819");
points(x=(cn4plus\$V9+cn4plus\$V10),y=(cn4plus\$V11),type="p",pch=19,cex=0.4,col="#E52420");
#label x
points(x=(cn2xchr\$V9+cn2xchr\$V10),y=(cn2xchr\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
points(x=(cn1xchr\$V9+cn1xchr\$V10),y=(cn1xchr\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
points(x=(cn3xchr\$V9+cn3xchr\$V10),y=(cn3xchr\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
points(x=(cn4xchr\$V9+cn4xchr\$V10),y=(cn4xchr\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
par(new=TRUE);
plot.default(x=c(1:10),y=c(1:10),ylim=c(0,110),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000");
lines(c(25,25),c(10,100),lty=2,col="black");
lines(((finalfactor * den2factor)),den2\$x,col="#67B32EAA",lwd=2);
lines(((finalfactor * den1factor)),den1\$x,col="#1C3660AA",lwd=2);
lines(((finalfactor * den3factor)),den3\$x,col="#F49819AA",lwd=2);
lines(((finalfactor * den4factor)),den4\$x,col="#E52420AA",lwd=2);
text(x=(finalfactor * cn1peakheight)+2,y=cn1peakpos,labels=signif(cn1peakpos,3),cex=0.7,srt=-90,col="#1C3660AA");
text(x=(finalfactor * cn3peakheight)+2,y=cn3peakpos,labels=signif(cn3peakpos,3),cex=0.7,srt=-90,col="#F49819AA");
text(x=(finalfactor * cn4peakheight)+2,y=cn4peakpos,labels=signif(cn4peakpos,3),cex=0.7,srt=-90,col="#E52420AA");
text(x=(finalfactor * cn2peakheight)+2,y=cn2peakpos,labels=signif(cn2peakpos,3),cex=0.7,srt=-90,col="#67B32EAA");
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(0,25),labels=c(0,sprintf("%.3f", maxden)),col="black",tck=0.01);
mtext("         CN Density",adj=0, cex=0.7, padj=-0.5);
par(mgp = c(3,1,0));
legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);

_END_OF_R_
        print R_COMMANDS "$R_command\n";
my $BRC = 0;
if ($BRC) {
        print R_COMMANDS "devoff <- dev.off();\n";

        #open up image for plotting
        if ($coverage_output_file =~ /.pdf/) {
            print R_COMMANDS "pdf(file=\"$coverage_output_file\",width=10,height=7.5,bg=\"white\");"."\n";
        }
        elsif ($coverage_output_file =~ /.png/) {
            print R_COMMANDS "png(file=\"$coverage_output_file\",width=1200,height=800);"."\n";
        }
        else {
            die "unrecognized coverage output file type...please append .pdf or .png to the end of your coverage output file\n";
        }

        print R_COMMANDS "par(mfrow=c(2,3));\n";
        print R_COMMANDS "absmaxx = 5000;\n";
}
        $R_command = <<"_END_OF_R_";

#transparency changes
#plotting coverage and copynumber on log scale with density plot incorporated with x chr labels
finalfactor = 25 / maxden100;
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
points(x=(cn2100x\$V9+cn2100x\$V10),y=(cn2100x\$V11),type="p",pch=19,cex=0.4,col="#67B32E");
points(x=(cn1minus100x\$V9+cn1minus100x\$V10),y=(cn1minus100x\$V11),type="p",pch=19,cex=0.4,col="#1C3660");
points(x=(cn3100x\$V9+cn3100x\$V10),y=(cn3100x\$V11),type="p",pch=19,cex=0.4,col="#F49819");
points(x=(cn4plus100x\$V9+cn4plus100x\$V10),y=(cn4plus100x\$V11),type="p",pch=19,cex=0.4,col="#E52420");
#label x
points(x=(cn2xchr100\$V9+cn2xchr100\$V10),y=(cn2xchr100\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
points(x=(cn1xchr100\$V9+cn1xchr100\$V10),y=(cn1xchr100\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
points(x=(cn3xchr100\$V9+cn3xchr100\$V10),y=(cn3xchr100\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
points(x=(cn4xchr100\$V9+cn4xchr100\$V10),y=(cn4xchr100\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
par(new=TRUE);
plot.default(x=c(1:10),y=c(1:10),ylim=c(0,110),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000");
lines(c(25,25),c(10,100),lty=2,col="black");
lines(((finalfactor * den2factor100)),den2100x\$x,col="#67B32EAA",lwd=2);
lines(((finalfactor * den1factor100)),den1100x\$x,col="#1C3660AA",lwd=2);
lines(((finalfactor * den3factor100)),den3100x\$x,col="#F49819AA",lwd=2);
lines(((finalfactor * den4factor100)),den4100x\$x,col="#E52420AA",lwd=2);
text(x=(finalfactor * cn1peakheight100)+2,y=cn1peakpos100,labels=signif(cn1peakpos100,3),cex=0.7,srt=-90,col="#1C3660AA");
text(x=(finalfactor * cn3peakheight100)+2,y=cn3peakpos100,labels=signif(cn3peakpos100,3),cex=0.7,srt=-90,col="#F49819AA");
text(x=(finalfactor * cn4peakheight100)+2,y=cn4peakpos100,labels=signif(cn4peakpos100,3),cex=0.7,srt=-90,col="#E52420AA");
text(x=(finalfactor * cn2peakheight100)+2,y=cn2peakpos100,labels=signif(cn2peakpos100,3),cex=0.7,srt=-90,col="#67B32EAA");
par(mgp = c(0, -1.4, 0));
axis(side=3,at=c(0,25),labels=c(0,sprintf("%.3f", maxden100)),col="black",tck=0.01);
par(mgp = c(3,1,0));
legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);
mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);

plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
points(y=additional_plot_points_cn1\$V2,x=additional_plot_points_cn1\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
points(y=additional_plot_points_cn2\$V2,x=additional_plot_points_cn2\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
points(y=additional_plot_points_cn3\$V2,x=additional_plot_points_cn3\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
plot.default(x=(z1\$V9+z1\$V10),y=(z1\$V11),log="x",xlab="Tumor Coverage",ylab="Tumor Variant Allele Frequency", main=paste(genome," Coverage"), type="p",pch=19,cex=0.4,col="#00000000",ylim=c(0,110),xlim=c(1,absmaxx));
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
points(y=additional_plot_points_cn4\$V2,x=additional_plot_points_cn4\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}




#############################################################


#final figure format
par(mfcol=c(5,3),mar=c(0.5,3,1,0),oma=c(3,0,4,0),mgp = c(3,1,0));
finalfactor = 25 / maxden100;

plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");
rect(0, 0, 100, 28, col = "#00000011",border=NA); #plot bg color
#lines(c(10,100),c(25,25),lty=2,col="black");
axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxden100)),las=1,cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);
lines(den2100x\$x,(finalfactor * den2factor100),col="#67B32EAA",lwd=2);
lines(den1100x\$x,(finalfactor * den1factor100),col="#1C3660AA",lwd=2);
lines(den3100x\$x,(finalfactor * den3factor100),col="#F49819AA",lwd=2);
lines(den4100x\$x,(finalfactor * den4factor100),col="#E52420AA",lwd=2);
text(x=cn1peakpos100,y=(finalfactor * cn1peakheight100)+1.7,labels=signif(cn1peakpos100,3),cex=0.7,srt=0,col="#1C3660AA");
text(x=cn3peakpos100,y=(finalfactor * cn3peakheight100)+1.7,labels=signif(cn3peakpos100,3),cex=0.7,srt=0,col="#F49819AA");
text(x=cn4peakpos100,y=(finalfactor * cn4peakheight100)+1.7,labels=signif(cn4peakpos100,3),cex=0.7,srt=0,col="#E52420AA");
text(x=cn2peakpos100,y=(finalfactor * cn2peakheight100)+1.7,labels=signif(cn2peakpos100,3),cex=0.7,srt=0,col="#67B32EAA");

axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
mtext(genome,adj=0,padj=-3.2,cex=0.65,side=3);


#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#00000055"); #plot bg color
#mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);
#legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);


#cn1plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(95,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn1minus100x\$V9+cn1minus100x\$V10),x=(cn1minus100x\$V11),type="p",pch=19,cex=0.4,col="#1C366044");
points(y=(cn1xchr100\$V9+cn1xchr100\$V10),x=(cn1xchr100\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn1\$V2,y=additional_plot_points_cn1\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 95, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#1C3660FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(1), cex=1, col="#FFFFFFFF") 

#cn2plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(95,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn2100x\$V9+cn2100x\$V10),x=(cn2100x\$V11),type="p",pch=19,cex=0.4,col="#67B32E44");
points(y=(cn2xchr100\$V9+cn2xchr100\$V10),x=(cn2xchr100\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn2\$V2,y=additional_plot_points_cn2\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 95, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#67B32EFF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(2), cex=1, col="#FFFFFFFF") 

#cn3plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(95,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn3100x\$V9+cn3100x\$V10),x=(cn3100x\$V11),type="p",pch=19,cex=0.4,col="#F4981999");
points(y=(cn3xchr100\$V9+cn3xchr100\$V10),x=(cn3xchr100\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn3\$V2,y=additional_plot_points_cn3\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 95, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#F49819FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(3), cex=1, col="#FFFFFFFF") 

#cn4plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(95,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn4plus100x\$V9+cn4plus100x\$V10),x=(cn4plus100x\$V11),type="p",pch=19,cex=0.4,col="#E5242044");
points(y=(cn4xchr100\$V9+cn4xchr100\$V10),x=(cn4xchr100\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn4\$V2,y=additional_plot_points_cn4\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 95, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#E52420FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(4), cex=1, col="#FFFFFFFF") 

axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.2);
mtext("Tumor Variant Allele Frequency",adj=0.5,padj=3.2,cex=0.5,side=1);



#all coverage points plotted
finalfactor = 25 / maxden;

plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");
rect(0, 0, 100, 28, col = "#00000011",border=NA); #plot bg color
#lines(c(10,100),c(25,25),lty=2,col="black");
axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxden)),las=1,cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);
lines(den2\$x,(finalfactor * den2factor),col="#67B32EAA",lwd=2);
lines(den1\$x,(finalfactor * den1factor),col="#1C3660AA",lwd=2);
lines(den3\$x,(finalfactor * den3factor),col="#F49819AA",lwd=2);
lines(den4\$x,(finalfactor * den4factor),col="#E52420AA",lwd=2);
text(x=cn1peakpos,y=(finalfactor * cn1peakheight)+1.7,labels=signif(cn1peakpos,3),cex=0.7,srt=0,col="#1C3660AA");
text(x=cn3peakpos,y=(finalfactor * cn3peakheight)+1.7,labels=signif(cn3peakpos,3),cex=0.7,srt=0,col="#F49819AA");
text(x=cn4peakpos,y=(finalfactor * cn4peakheight)+1.7,labels=signif(cn4peakpos,3),cex=0.7,srt=0,col="#E52420AA");
text(x=cn2peakpos,y=(finalfactor * cn2peakheight)+1.7,labels=signif(cn2peakpos,3),cex=0.7,srt=0,col="#67B32EAA");

axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
mtext(genome,adj=0,padj=-3.2,cex=0.65,side=3);


#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#00000055"); #plot bg color
#mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);
#legend(x="topright",horiz=TRUE,xjust=0, c("1", "2", "3", "4+","Chr $chr_highlight"),col=c("#1C3660","#67B32E","#F49819","#E52420","#A020F0"),pch=c(19,19,19,19,2),cex=0.6);


#cn1plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(5,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn1minus\$V9+cn1minus\$V10),x=(cn1minus\$V11),type="p",pch=19,cex=0.4,col="#1C366044");
points(y=(cn1xchr\$V9+cn1xchr\$V10),x=(cn1xchr\$V11),type="p",pch=2,cex=0.8,col="#1C366044");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn1\$V2,y=additional_plot_points_cn1\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#1C3660FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(1), cex=1, col="#FFFFFFFF") 

#cn2plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(5,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn2\$V9+cn2\$V10),x=(cn2\$V11),type="p",pch=19,cex=0.4,col="#67B32E44");
points(y=(cn2xchr\$V9+cn2xchr\$V10),x=(cn2xchr\$V11),type="p",pch=2,cex=0.8,col="#67B32E44");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn2\$V2,y=additional_plot_points_cn2\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#67B32EFF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(2), cex=1, col="#FFFFFFFF") 

#cn3plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(5,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn3\$V9+cn3\$V10),x=(cn3\$V11),type="p",pch=19,cex=0.4,col="#F4981999");
points(y=(cn3xchr\$V9+cn3xchr\$V10),x=(cn3xchr\$V11),type="p",pch=2,cex=0.8,col="#F4981955");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn3\$V2,y=additional_plot_points_cn3\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#F49819FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(3), cex=1, col="#FFFFFFFF") 

#cn4plot
plot.default(x=(z1\$V11),y=(z1\$V9+z1\$V10),log="y", type="p",pch=19,cex=0.4,col="#00000000",xlim=c(-1,101),ylim=c(5,absmaxx+5),axes=FALSE, ann=FALSE,xaxs="i",yaxs="i");
points(y=(cn4plus\$V9+cn4plus\$V10),x=(cn4plus\$V11),type="p",pch=19,cex=0.4,col="#E5242044");
points(y=(cn4xchr\$V9+cn4xchr\$V10),x=(cn4xchr\$V11),type="p",pch=2,cex=0.8,col="#E5242044");
#add in highlight of points selected for by script input
if(length(additional_plot_points) > 1) {
	points(x=additional_plot_points_cn4\$V2,y=additional_plot_points_cn4\$V3,type="p",pch=7,cex=0.8,col="#555555FF");
}
axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
for (i in 2:length(axTicks(2)-1)) {
	lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
}
rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA); #plot bg color
#add cn circle
points(x=c(97),y=c((absmaxx+5)*0.70),type="p",pch=19,cex=3,col="#E52420FF");
text(c(97),y=c((absmaxx+5)*0.70), labels=c(4), cex=1, col="#FFFFFFFF") 

axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.2);
mtext("Tumor Variant Allele Frequency",adj=0.5,padj=3.2,cex=0.5,side=1);

devoff <- dev.off();








#copy number pdf starts here
pdf(file=\"$copynumber_output_file\",width=10,height=7.5);
par(mfrow=c(2,2));
cov1=(z1\$V20);
cov2=(z2\$V20);
maxx3=max(c(cov1,cov2));
if (maxx <= 4) {maxx = 4};

plot.default(((den2\$y * 100)+30),den2\$x,xlab="Normal Variant Allele Frequency",ylab="Tumor Variant Allele Frequency",  main=paste(genome," Variant Allele Frequency"),xlim=c(0,100),ylim=c(0,100),col="#67B32E44", type="p",pch=19,cex=0.4);
points(((den1\$y * 100)+20),den1\$x,col="#E5242055", type="p",pch=19,cex=0.4);
points(((den3\$y * 100)+40),den3\$x,col="#1C366044", type="p",pch=19,cex=0.4);
points(((den4\$y * 100)+50),den4\$x,col="#F49819FF", type="p",pch=19,cex=0.4);
points(x=cn1minus\$V7,y=cn1minus\$V11, type="p",pch=19,cex=0.4,col="#E5242055");
points(x=cn2\$V7,y=cn2\$V11, type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=cn3\$V7,y=cn3\$V11, type="p",pch=19,cex=0.4,col="#1C366044");
points(x=cn4plus\$V7,y=cn4plus\$V11, type="p",pch=19,cex=0.4,col="#F49819FF");
legend(x="topright", title = "Copy Number", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),pch=19);
legend(x="bottomright", title = "Copy Number Density", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),lty = c(1,1,1), lwd = 1);

plot.default(x=cn1minus100x\$V7,y=cn1minus100x\$V11,xlab="Normal Variant Allele Frequency",ylab="Tumor Variant Allele Frequency", main=paste(genome," Variant Allele Frequency"), type="p",pch=19,cex=0.4,col="#E5242055",xlim=c(0,100),ylim=c(0,100));
points(x=cn2100x\$V7,y=cn2100x\$V11, type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=cn3100x\$V7,y=cn3100x\$V11, type="p",pch=19,cex=0.4,col="#1C366044");
points(x=cn4plus100x\$V7,y=cn4plus100x\$V11, type="p",pch=19,cex=0.4,col="#F49819FF");
points(((den1100x\$y * 100)+20),den1100x\$x,col="#E5242055", type="p",pch=19,cex=0.4);
points(((den2100x\$y * 100)+30),den2100x\$x,col="#67B32E44", type="p",pch=19,cex=0.4);
points(((den3100x\$y * 100)+40),den3100x\$x,col="#1C366044", type="p",pch=19,cex=0.4);
points(((den4100x\$y * 100)+50),den4100x\$x,col="#F49819FF", type="p",pch=19,cex=0.4);
legend(x="topright", title = "Copy Number", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),pch=19);
legend(x="bottomright", title = "Copy Number Density", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),lty = c(1,1,1), lwd = 1);
mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);

plot.default(x=((cn1minus\$V7*cn1minus\$V21)-100),y=((cn1minus\$V11*cn1minus\$V20)-100),xlab="Normal Variant Allele Frequency (CN Corrected)",ylab="Tumor Variant Allele Frequency (CN Corrected)", main=paste(genome," Copy Number"), type="p",pch=19,cex=0.4,col="#E5242055",xlim=c(-100,100),ylim=c(-100,100));
points(x=((cn2\$V7*cn2\$V21) - 100),y=((cn2\$V11*cn2\$V20)-100), type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=((cn3\$V7*cn3\$V21)-100),y=((cn3\$V11*cn3\$V20)-100), type="p",pch=19,cex=0.4,col="#1C366044");
points(x=((cn4plus\$V7*cn4plus\$V21)-100),y=((cn4plus\$V11*cn4plus\$V20)-100), type="p",pch=19,cex=0.4,col="#F49819FF");
legend(x="topright", title = "Copy Number", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),pch=19);

plot.default(x=((cn1minus100x\$V7*cn1minus100x\$V21)-100),y=((cn1minus100x\$V11*cn1minus100x\$V20)-100),xlab="Normal Variant Allele Frequency (CN Corrected)",ylab="Tumor Variant Allele Frequency (CN Corrected)", main=paste(genome," Copy Number"), type="p",pch=19,cex=0.4,col="#E5242055",xlim=c(-100,100),ylim=c(-100,100));
points(x=((cn2100x\$V7*cn2100x\$V21) - 100),y=((cn2100x\$V11*cn2100x\$V20)-100), type="p",pch=19,cex=0.4,col="#67B32E44");
points(x=((cn3100x\$V7*cn3100x\$V21)-100),y=((cn3100x\$V11*cn3100x\$V20)-100), type="p",pch=19,cex=0.4,col="#1C366044");
points(x=((cn4plus100x\$V7*cn4plus100x\$V21)-100),y=((cn4plus100x\$V11*cn4plus100x\$V20)-100), type="p",pch=19,cex=0.4,col="#F49819FF");
legend(x="topright", title = "Copy Number", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),pch=19);
mtext("Normal and Tumor Coverage > 100",cex=0.7, padj=-0.5);

par(mfrow=c(1,1));
s3d <- scatterplot3d(x=cov100xplus\$V7,z=cov100xplus\$V11,y=cov100xplus\$V20, type="p", angle=55, scale.y=0.7, cex.symbols=0.4, pch=19,xlab="Normal Variant Allele Frequency",zlab="Tumor Variant Allele Frequency",ylab="Copy Number",xlim=c(0,100),zlim=c(0,100),ylim=c(0,5),color="#67B32EFF",box=FALSE);
s3d\$points3d(x=cov20x\$V7,z=cov20x\$V11,y=cov20x\$V20, type="p",pch=19,cex=0.4,col="#E52420FF");
s3d\$points3d(x=cov50x\$V7,z=cov50x\$V11,y=cov50x\$V20, type="p",pch=19,cex=0.4,col="#1C3660FF");
s3d\$points3d(x=cov100x\$V7,z=cov100x\$V11,y=cov100x\$V20, type="p",pch=19,cex=0.4,col="#F49819FF");
legend(x="topright", title = "Coverage", c("0-20x", "20-50x", "50-100x", "100+x"),col=c("#E52420","#1C3660","#F49819","#67B32E"),pch=19);

cn1cov = (cn1minus\$V9+cn1minus\$V10);
cn2cov = (cn2\$V9+cn2\$V10);
cn3cov = (cn3\$V9+cn3\$V10);
cn4cov = (cn4plus\$V9+cn4plus\$V10);
s3d <- scatterplot3d(x=cn2\$V7,z=cn2\$V11,y=cn2cov, type="p", angle=55, scale.y=0.7, cex.symbols=0.4, pch=19,xlab="Normal Variant Allele Frequency",zlab="Tumor Variant Allele Frequency",ylab="Coverage",xlim=c(0,100),zlim=c(0,100),ylim=c(0,maxx),color="#67B32EFF",box=FALSE);
s3d\$points3d(x=cn1minus\$V7,z=cn1minus\$V11,y=cn1cov, type="p",pch=19,cex=0.4,col="#E52420FF");
s3d\$points3d(x=cn3\$V7,z=cn3\$V11,y=cn3cov, type="p",pch=19,cex=0.4,col="#1C3660FF");
s3d\$points3d(x=cn4plus\$V7,z=cn4plus\$V11,y=cn4cov, type="p",pch=19,cex=0.4,col="#F49819FF");
legend(x="topright", title = "Copy Number", c("1", "2", "3", "4+"),col=c("#E52420","#67B32E","#1C3660","#F49819"),pch=19);

#plot.default(den,xlab="Normal Variant Allele Frequency",ylab="Tumor Variant Allele Frequency",  main=paste(genome," Variant Allele Frequency"));
#library(lattice);
#cloud(cn1minus\$V20 ~ cn1minus\$V7 * cn1minus\$V11);
#cloud(cn2\$V20 ~ cn2\$V7 * cn2\$V11);
#cloud(cn3\$V20 ~ cn2\$V7 * cn3\$V11);
#cloud(cn4plus\$V20 ~ cn4plus\$V7 * cn4plus\$V11);

devoff <- dev.off();

q();

_END_OF_R_
        print R_COMMANDS "$R_command\n";

	close R_COMMANDS;

	my $cmd = "R --vanilla --slave \< $r_script_output_file";
	my $return = Genome::Sys->shellcmd(
           cmd => "$cmd",
           output_files => [$coverage_output_file, $copynumber_output_file],
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
    my ($file, $tumnor)=@_;
    my %info_hash;
    my $fh=new FileHandle($file);
    while(my $line = <$fh>)
    {
	chomp($line);
	unless ($line =~ /^\w+\t\d+\t\d+\t/) { next;}
	my ($chr,$start,$end,$size,$nmarkers,$cn,$adjusted_cn,$cn_normal,$adjusted_cn_normal,$score)=split(/\t/, $line);
	my $pos=$start."_".$end;
	if ($tumnor eq 'tumor') {
		$info_hash{$chr}{$pos}=$adjusted_cn;
	}
	elsif ($tumnor eq 'normal') {
		$info_hash{$chr}{$pos}=$adjusted_cn_normal;
	}
    }
    return \%info_hash;
}




