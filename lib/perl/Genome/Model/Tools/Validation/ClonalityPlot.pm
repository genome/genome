package Genome::Model::Tools::Validation::ClonalityPlot;

use strict;
use warnings;
use FileHandle;
use Genome;
use File::Basename;
use FileHandle;

class Genome::Model::Tools::Validation::ClonalityPlot {
    is => 'Command',                       

    has => [
    varscan_file => { 
        is => 'Text', 
        doc => "File of varscan validated calls, ex: ", 
        is_optional => 0, 
        is_input => 1 },

    varscan_file_with_cn => {
        is => 'Text',
        doc => "Optional output file of varscan validated calls with copy-number appended",
        is_optional => 1 },

    cnvhmm_file => { 
        is => 'Text',
        doc => "File of cnvhmm whole genome predictions",
        is_optional => 1,
        is_input => 1 },        

    cbs_file => { 
        is => 'Text',
        doc => "File of CN-altered segments according to CBS - assumes segment mean is log2 value",
        is_optional => 1,
        is_input => 1 },        

    sample_id => {
        is => 'Text',
        doc => "Sample ID to be put on graphs",
        is_optional => 1,
        is_input => 1,
        default => 'unspecified' },

    analysis_type => {
        is => 'Text',
        doc => "Either \'wgs\' for somatic pipeline output or \'capture\' for validation pipeline output",
        is_optional => 1,
        is_input => 1,
        default => 'capture'},

    chr_highlight => {
        is => 'Text',
        doc => "Choose a Chromosome to Highlight with Purple Circles on Plot",
        is_optional => 1,
        is_input => 1,
        default => 'X'},

    positions_highlight => {
        is => 'Text',
        doc => "A tab-delim file list of positions chr\\tposition to highlight on plots",
        is_optional => 1,
        is_input => 1},

    r_script_output_file => {
        is => 'Text',
        doc => "R script built and run by this module",
        is_optional => 0,
        is_input => 1},

    output_image => {
        is => 'Text',
        doc => "PDF Coverage output file",
        is_optional => 0,
        is_input => 1,
        is_output => 1 },

    skip_if_output_is_present => {
        is => 'Text',
        doc => "Skip if Output is Present",
        is_optional => 1,
        is_input => 1,
        default => 0},

    minimum_labelled_peak_height => {
        is => 'Text',
        doc => "only peaks that exceed this height get labelled",
        is_optional => 1,
        is_input => 1,
        default => 0.001},        

    only_label_highest_peak => {
        is => 'Boolean',
        doc => "only label the highest peak",
        is_optional => 1,
        default => 0},        

    plot_only_cn2 => {
        is => 'Boolean',
        doc => "only plot the CN2 data",
        is_optional => 1,
        default => 0},        

    plot_clusters => {
        is => 'Boolean',
        doc => "plot clustered points on the CN2 scatterplot",
        is_optional => 1,
        default => 0},

    tumor_purity => {
        is => 'Integer',
        doc => "tumor purity in terms of a percent from 1 to 100; will be estimated by the tool as a default, but this may not be super-reliable",
        is_optional => 1,
        default => 0},

    clustered_data_output_file => {
        is => 'String',
        doc => "If clustering data, this output file will contain the points used in clustering and their associated clusters",
        is_optional => 1,
        default => "NULL"},

    max_clusters_to_look_for => {
        is => 'Integer',
        doc => "If clustering data, this is the maximum number of subclones to fit components to within the data",
        is_optional => 1,
        default => 8},

    component_distribution => {
        is => 'String',
        doc => "If clustering data, this will be the distribution used for fitting components. Please choose between 'Binomial' or 'Normal'.",
        is_optional => 1,
        default => "Binomial"},
    ],
};

sub sub_command_sort_position { 1 }

sub help_brief {
    "Plot CN-separated SNV density plots with optional clustering"
}

sub help_synopsis {
    return <<EOS
        Inputs of Varscan and copy-number segmentation data, Output of R plots.
EXAMPLE:	gmt validation clonality-plot --analysis-type 'capture' --varscan-file snvs.txt --cnvhmm-file cnaseg.txt --output-image clonality.pdf --sample-id 'Sample'
EXAMPLE:	gmt validation clonality-plot --analysis-type 'capture' --varscan-file snvs.txt --cnvhmm-file cnaseg.txt --output-image clonality.pdf --sample-id 'Sample' --r-script-output-file Sample.R --positions-highlight chr.pos.txt --varscan-file-with-cn snvs.txt.cn
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS
This tool can be used to plot the CN-separated SNV density plots that are known at TGI as 'clonality plots'. Can be used for WGS or Capture data, but was mostly intended for Capture data, and hence the SNV file format is currently Varscan output.
EOS
}

sub execute {
    my $self = shift;

    ##inputs##
    my $varscan_file = $self->varscan_file;
    my $varscan_file_with_cn = $self->varscan_file_with_cn;
    my $cnvhmm_file = $self->cnvhmm_file;
    my $cbs_file = $self->cbs_file;
    my $sample_id = $self->sample_id;
    my $readcount_cutoff;
    my $chr_highlight = $self->chr_highlight;
    my $positions_highlight = $self->positions_highlight;
    my $tumor_purity = $self->tumor_purity;
    my $component_dist = $self->component_distribution;
    my $max_clusters_to_look_for = $self->max_clusters_to_look_for;
    ##outputs##
    my $r_script_output_file = $self->r_script_output_file;
    my $output_image = $self->output_image;
    my $clustered_data_output_file = $self->clustered_data_output_file;
    ##options##
    my $skip_if_output_is_present = $self->skip_if_output_is_present;
    my $analysis_type = $self->analysis_type;
    my $minimum_labelled_peak_height = $self->minimum_labelled_peak_height;
    my $only_label_highest_peak = $self->only_label_highest_peak;
    my $plot_only_CN2 = $self->plot_only_cn2;
    my $plot_clusters = $self->plot_clusters;

    #set readcount cutoffs
    if ($analysis_type eq 'wgs') {
        $readcount_cutoff = 20;
    }
    elsif ($analysis_type eq 'capture') {
        $readcount_cutoff = 100;
    }
    else {
        die "analysis type: $analysis_type not supported, choose either wgs or capture";
    }


    #create a hash of positions to be highlighted
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
    my ($tfh,$temp_path);
    if (defined $self->varscan_file_with_cn) {
        $temp_path = $self->varscan_file_with_cn;
        $tfh = FileHandle->new($temp_path,"w");
    }
    else {
        ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    }

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

    my %copynumber_hash_tumor;
    if(defined($cnvhmm_file)){
        #build the copy number hashes
        %copynumber_hash_tumor=%{&build_hash($cnvhmm_file)};
    } elsif(defined($cbs_file)){
        #build the copy number hashes
        %copynumber_hash_tumor=%{&build_hash_cbs($cbs_file)};
    } 


    #read in the varscan input
    my $varscan_input = new FileHandle ($varscan_file);
    while (my $line2 = <$varscan_input>) {
        chomp($line2);
        my ($chr, $pos, $ref, $var, 
            $normal_ref, $normal_var, $normal_var_pct, $normal_IUB, 
            $tumor_ref, $tumor_var, $tumor_var_pct, $tumor_IUB, 
            $varscan_call, $germline_pvalue, $somatic_pvalue, @otherstuff) = split(/\t/, $line2);

        my $varscan_cn_tumor=&get_cn($chr,$pos,$pos,\%copynumber_hash_tumor);

        print $tfh "$line2\t$varscan_cn_tumor\n";
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

    my $dir_name = dirname(__FILE__);
    my $r_library = $dir_name . "/VarScanGraphLib.R";

    #create the R script that will be run to produce the plot
#-------------------------------------------------
    my $R_command = <<"_END_OF_R_";
#options(echo = FALSE);#suppress output to stdout
#sink("/dev/null");
    genome=paste(\"$sample_id\","Clonality Plot",sep=\" \");
    source(\"$r_library\"); #this contains R functions for loading and graphing VarScan files

    library(fpc);
    library(scatterplot3d);
    varscan.load_snp_output(\"$temp_path\",header=F)->xcopy;
    varscan.load_snp_output(\"$temp_path\",header=F,min_tumor_depth=$readcount_cutoff,min_normal_depth=$readcount_cutoff)->xcopy100;
    clustered_data_output_file = \"$clustered_data_output_file\";

    additional_plot_points = 0;
    additional_plot_points_cn1 = 0;
    additional_plot_points_cn2 = 0;
    additional_plot_points_cn3 = 0;
    additional_plot_points_cn4 = 0;
_END_OF_R_
#-------------------------------------------------

    print R_COMMANDS "source(\"" . $dir_name . "/ClonalityPlot.R\")\n";


    print R_COMMANDS "$R_command\n";


    #add highlighting code if specified
    if ($positions_highlight && -s $positions_highlight && 1 && 1) { #these &&1 mean nothing, they just make my text editor color things correctly (it hates -s without being s///)
#-------------------------------------------------
        $R_command = <<"_END_OF_R_";
        additional_plot_points <- read.table(\"$temp_path2\", header = FALSE, sep = "\t");
        additional_plot_points_cn1=subset(additional_plot_points, additional_plot_points\$V1 >= 0 & additional_plot_points\$V1 <= 1.75);
        additional_plot_points_cn2=subset(additional_plot_points, additional_plot_points\$V1 >= 1.75 & additional_plot_points\$V1 <= 2.25);
        additional_plot_points_cn3=subset(additional_plot_points, additional_plot_points\$V1 >= 2.25 & additional_plot_points\$V1 <= 3.5);
        additional_plot_points_cn4=subset(additional_plot_points, additional_plot_points\$V1 >= 3.5);
_END_OF_R_
#-------------------------------------------------

        print R_COMMANDS "$R_command\n";
    }


#-------------------------------------------------
    $R_command = <<"_END_OF_R_";    
    z1=subset(xcopy, xcopy\$V13 == "Somatic");
    z2=subset(xcopy100, xcopy100\$V13 == "Somatic");
    xchr=subset(z1,z1\$V1 == "$chr_highlight");
    xchr100=subset(z2,z2\$V1 == "$chr_highlight");
    covtum1=(z1\$V9+z1\$V10);
    covtum2=(z2\$V9+z2\$V10);
    absmaxx=maxx=max(c(covtum1,covtum2));
#    covnorm1=(z1\$V5+z1\$V6);
#    covnorm2=(z2\$V5+z2\$V6);
#    absmaxx2=maxx2=max(c(covnorm1,covnorm2));

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

    if(dim(cn1minus)[1] < 2) {
        den1\$x = den1\$y=1000;
    } else {
        den1 <- density(cn1minus\$V11, from=0,to=100,na.rm=TRUE); den1factor = dim(cn1minus)[1]/N * den1\$y;}
    if(dim(cn2)[1] < 2) {
        den2\$x = den2\$y=1000;
    } else {
        den2 <- density(cn2\$V11, from=0,to=100,na.rm=TRUE); den2factor = dim(cn2)[1]/N * den2\$y;};
    if(dim(cn3)[1] < 2) {
        den3\$x = den3\$y=1000;
    } else {
        den3 <- density(cn3\$V11, from=0,to=100,na.rm=TRUE); den3factor = dim(cn3)[1]/N * den3\$y;};
    if(dim(cn4plus)[1] < 2) {
        den4\$x = den4\$y=1000;
    } else {
        den4 <- density(cn4plus\$V11, from=0,to=100,na.rm=TRUE); den4factor = dim(cn4plus)[1]/N * den4\$y;};
    if(dim(cn1minus100x)[1] < 2) {
        den1100x\$x = den1100x\$y=1000
    } else {
        den1100x <- density(cn1minus100x\$V11, from=0,to=100,na.rm=TRUE); den1factor100 = dim(cn1minus100x)[1]/N100 * den1100x\$y;};
    if(dim(cn2100x)[1] < 2) {
        den2100x\$x = den2100x\$y=1000
    } else {
        den2100x <- density(cn2100x\$V11, from=0,to=100,na.rm=TRUE);den2factor100 = dim(cn2100x)[1]/N100 * den2100x\$y;};
    if(dim(cn3100x)[1] < 2) {
        den3100x\$x = den3100x\$y=1000
    } else {
        den3100x <- density(cn3100x\$V11, from=0,to=100,na.rm=TRUE);den3factor100 = dim(cn3100x)[1]/N100 * den3100x\$y;};
    if(dim(cn4plus100x)[1] < 2) {
        den4100x\$x = den4100x\$y=1000
    } else {
        den4100x <- density(cn4plus100x\$V11, from=0,to=100,na.rm=TRUE);den4factor100 = dim(cn4plus100x)[1]/N100 * den4100x\$y
    }

    # dennormcov <- density((z1\$V5+z1\$V6), bw=4, from=0,to=maxx,na.rm=TRUE);
    # dennormcov100x <- density((z2\$V5+z2\$V6), bw=4, from=0,to=maxx,na.rm=TRUE);
    # dentumcov <- density((z1\$V9+z1\$V10), bw=4, from=0,to=maxx,na.rm=TRUE);
    # dentumcov100x <- density((z2\$V9+z2\$V10), bw=4, from=0,to=maxx,na.rm=TRUE);

#find inflection points (peaks)
#den2diff = diff(den2\$y);

    peaks<-function(series,span=3){
        z <- embed(series, span);
        s <- span%/%2;
        v<- max.col(z) == 1 + s;
        result <- c(rep(FALSE,s),v);
        result <- result[1:(length(result)-s)];
        result;
    } 

    #labels to use for density plot values
    if(dim(cn1minus)[1] < 2) {
        cn1peaks = cn1peakpos = cn1peakheight = 0;
    } else {
        cn1peaks = peaks(den1factor); 
        cn1peaks = append(cn1peaks,c("FALSE","FALSE"),after=length(cn1peaks)); 
        cn1peakpos = subset(den1\$x,cn1peaks==TRUE & den1\$y > 0.001); 
        cn1peakheight = subset(den1factor,cn1peaks==TRUE & den1\$y > 0.001);
    }

    if(dim(cn2)[1] < 2) {
        cn2peaks = cn2peakpos = cn2peakheight = 0;
    } else {
        cn2peaks = peaks(den2factor);
        cn2peaks = append(cn2peaks,c("FALSE","FALSE"),after=length(cn2peaks)); 
        cn2peakpos = subset(den2\$x,cn2peaks==TRUE & den2\$y > 0.001); 
        ###print(cn2peakpos);
        cn2peakheight = subset(den2factor,cn2peaks==TRUE & den2\$y > 0.001);
    }

    if(dim(cn3)[1] < 2) {
        cn3peaks = cn3peakpos = cn3peakheight = 0;
    } else {
        cn3peaks = peaks(den3factor);
        cn3peaks = append(cn3peaks,c("FALSE","FALSE"),after=length(cn3peaks));
        cn3peakpos = subset(den3\$x,cn3peaks==TRUE & den3\$y > 0.001);
        cn3peakheight = subset(den3factor,cn3peaks==TRUE & den3\$y > 0.001);}
    if(dim(cn4plus)[1] < 2) {cn4peaks = cn4peakpos = cn4peakheight = 0;} else {cn4peaks = peaks(den4factor);
        cn4peaks = append(cn4peaks,c("FALSE","FALSE"),after=length(cn4peaks));
        cn4peakpos = subset(den4\$x,cn4peaks==TRUE & den4\$y > 0.001);
        cn4peakheight = subset(den4factor,cn4peaks==TRUE & den4\$y > 0.001);
    }
    if(dim(cn1minus100x)[1] < 2) {
        cn1peaks100 = cn1peakpos100 = cn1peakheight100 = 0;
    } else {
        cn1peaks100 = peaks(den1factor100);
        cn1peaks100 = append(cn1peaks100,c("FALSE","FALSE"),after=length(cn1peaks100));
        cn1peakpos100 = subset(den1100x\$x,cn1peaks100==TRUE & den1100x\$y > 0.001);
        cn1peakheight100 = subset(den1factor100,cn1peaks100==TRUE & den1100x\$y > 0.001);}
    if(dim(cn2100x)[1] < 2) {
        cn2peaks100 = cn2peakpos100 = cn2peakheight100 = 0;
    } else {
        cn2peaks100 = peaks(den2factor100);
        cn2peaks100 = append(cn2peaks100,c("FALSE","FALSE"),after=length(cn2peaks100));
        cn2peakpos100 = subset(den2100x\$x,cn2peaks100==TRUE & den2100x\$y > 0.001);
        cn2peakheight100 = subset(den2factor100,cn2peaks100==TRUE & den2100x\$y > 0.001);}
    if(dim(cn3100x)[1] < 2) {
        cn3peaks100 = cn3peakpos100 = cn3peakheight100 = 0;
    } else {
        cn3peaks100 = peaks(den3factor100);
        cn3peaks100 = append(cn3peaks100,c("FALSE","FALSE"),after=length(cn3peaks100));
        cn3peakpos100 = subset(den3100x\$x,cn3peaks100==TRUE & den3100x\$y > 0.001);
        cn3peakheight100 = subset(den3factor100,cn3peaks100==TRUE & den3100x\$y > 0.001);}
    if(dim(cn4plus100x)[1] < 2) {
        cn4peaks100 = cn4peakpos100 = cn4peakheight100 = 0;
    } else {
        cn4peaks100 = peaks(den4factor100);
        cn4peaks100 = append(cn4peaks100,c("FALSE","FALSE"),after=length(cn4peaks100));
        cn4peakpos100 = subset(den4100x\$x,cn4peaks100==TRUE & den4100x\$y > 0.001);
        cn4peakheight100 = subset(den4factor100,cn4peaks100==TRUE & den4100x\$y > 0.001);
    }
    maxden = max(c(den1factor,den2factor,den3factor,den4factor));
    maxden100 = max(c(den1factor100,den2factor100,den3factor100,den4factor100));

    #determine tumor purity
    analysis_type = \"$analysis_type\";
    purity = $tumor_purity;
    if (purity == 0) {
    if (analysis_type == \"wgs\") {
        purity = max(cn2peakpos[which(cn2peakpos<=50)])*2;
        next_highest_peak = min(cn2peakpos[which(cn2peakpos>50)]);
        if (next_highest_peak > 50 && next_highest_peak < 60 && purity < 60) { purity = 100; }
    } else { #capture
        purity = max(cn2peakpos100[which(cn2peakpos100<=50)])*2;
        next_highest_peak = min(cn2peakpos100[which(cn2peakpos100>50)]);
        if (next_highest_peak > 50 && next_highest_peak < 60 && purity < 60) { purity = 100; }
    }
    if (purity == 0) { purity = 100; }
    print(paste("Tumor purity estimated to be ",purity,".",sep=""));
    }

    #determine number of clusters in the dataset
    num_clusters = 0;
    if ($plot_clusters) {

    print("Performing 'mixdist' analysis...");
    max_clusters_to_look_for = $max_clusters_to_look_for;
    component_dist = \"$component_dist\";

    #function to process mixdist results
    process_percents <- function(percents,chisq,pval,distr) {
    minchisq = NULL;
    best_fit_component_assessment = 0;
    minpval = NULL;
    true_cluster_count = NULL;

    for (i in 1:max_clusters_to_look_for) {
    if (percents[i]!="Error" && !is.nan(pval[i]) && pval[i] < 0.05) {
    if (is.null(minchisq) || chisq[i] < minchisq) {
    minchisq = chisq[i];
    minpval = pval[i];
    best_fit_component_assessment = i;
    }
    }
    true_cluster_count[i] = 0;
    percentage_per_cluster = as.numeric(percents[[i]]);
    true_clusters = percentage_per_cluster > 0.02;
    for (j in 1:i) {
    if (isTRUE(true_clusters[j])) {
    true_cluster_count[i] = true_cluster_count[i] + 1;
    }
    }
    }
    print(paste("chosen_component_assessment = ",best_fit_component_assessment,", true_component_count = ",true_cluster_count[best_fit_component_assessment],", assessment_chisq_value = ",minchisq,", assessment_p-value = ",minpval,sep=""));
    return(true_cluster_count[best_fit_component_assessment]);
    }

    # run mixdist
    library(mixdist);
    data = NULL;
    if (analysis_type == \"wgs\") { data = round(cn2\$V11); } else { data = round(cn2100x\$V11); }
    grouped_data = NULL;

    # if ceiling of max(data) is odd, then add 1 and group data. else, group data using the even ceiling
    if (ceiling(max(data))%%2) {
    grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data))+1,2)));
    } else { 
    grouped_data = mixgroup(data, breaks = c(seq(0,ceiling(max(data)),2)));
    }

    # for each component count, get mixdist fit estimates for normal distribution
    percents=NULL; chisq=NULL; pval=NULL; distr=component_dist;
    for (i in 1:max_clusters_to_look_for) {

    data_params = NULL;
    test = NULL;

    if (distr == \"Normal\") {
    data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data),i));
    test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="norm"), silent=TRUE);
    }
    if (distr == \"Binomial\") {
    data_params = mixparam(c(1:i)*(purity/2)/i,rep(sd(data)/2,i));
    test=try(mix(mixdat=grouped_data,mixpar=data_params, emsteps=3, dist="binom", constr=mixconstr(consigma="BINOM",size=rep(round(length(data)),i))), silent=TRUE);
    }

    if (class(test) == 'try-error') {
    percents[i] = "Error";
    print(paste("Common mixdist error when looking for ",i," components.",sep=""));
    }
    else {
    percents[i]=list(test\$parameters\$pi)
    chisq[i]=test\$chisq;
    pval[i] = test\$P;

    #print(paste("PLOTTING ",i)); filename = paste("plot_component_",i,".pdf",sep=""); dev.new(); plot(test); dev.copy(pdf,filename); dev.off(); ## TEST
    }
    }
    num_clusters = process_percents(percents,chisq,pval,distr);
    }

_END_OF_R_
#-------------------------------------------------

print R_COMMANDS "$R_command\n";


#open up image for plotting
if ($output_image =~ /.pdf/) {
print R_COMMANDS "pdf(file=\"$output_image\",width=3.3,height=7.5,bg=\"white\");"."\n";

} elsif ($output_image =~ /.png/) {
    print R_COMMANDS "png(file=\"$output_image\",width=400,height=800);"."\n";

} else {
    die "unrecognized coverage output file type...please append .pdf or .png to the end of your coverage output file\n";
}
print R_COMMANDS "par(mfcol=c(5,1),mar=c(0.5,3,1,1.5),oma=c(3,0,4,0),mgp = c(3,1,0));"."\n";
# } else {
#     if ($output_image =~ /.pdf/) {
#         print R_COMMANDS "pdf(file=\"$output_image\",width=3.3,height=3.3,bg=\"white\");"."\n";

#     } elsif ($output_image =~ /.png/) {
#         print R_COMMANDS "png(file=\"$output_image\",width=400,height=340);"."\n";

#     } else {
#         die "unrecognized coverage output file type...please append .pdf or .png to the end of your coverage output file\n";
#     }

#     print R_COMMANDS "par(mfcol=c(2,1),mar=c(0.5,3,1,1.5),oma=c(3,0,4,0),mgp = c(3,1,0));"."\n";
# }


if ($analysis_type eq 'capture') {
#-------------------------------------------------
    $R_command = <<"_END_OF_R_";
    #final figure format
    finalfactor = 25 / maxden100;

    plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");
    rect(0, 0, 100, 28, col = "#00000011",border=NA); #plot bg color
    #lines(c(10,100),c(25,25),lty=2,col="black");
    axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxden100)),las=1,cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);
    lines(den2100x\$x,(finalfactor * den2factor100),col="#67B32EAA",lwd=2);
    if(!($plot_only_CN2)){
        lines(den1100x\$x,(finalfactor * den1factor100),col="#1C3660AA",lwd=2)\n;
        lines(den3100x\$x,(finalfactor * den3factor100),col="#F49819AA",lwd=2)\n;
        lines(den4100x\$x,(finalfactor * den4factor100),col="#E52420AA",lwd=2)\n;
    }

    #legend
    if ($plot_only_CN2) {
        legend(x="topright",lwd=2,legend=c("2 Copies"), col=c("#67B32EAA"),bty="n",cex=0.6, y.intersp=1.25);
    } else { 
        legend(x="topright",lwd=2,legend=c("1 Copy","2 Copies","3 Copies","4 Copies"), col=c("#1C3660AA","#67B32EAA","#F49819AA","#E52420AA"),bty="n",cex=0.6, y.intersp=1.25);
    }

    if(!($plot_only_CN2)){
        ppos = c();
        if($only_label_highest_peak){
            ppos = which((cn1peakheight100 == max(cn1peakheight100)) & (cn1peakheight100 > $minimum_labelled_peak_height));
        } else {
            ppos = which(cn1peakheight100 > $minimum_labelled_peak_height);
        }
        if(!(length(ppos) == 0)){
            text(x=cn1peakpos100[ppos],
                 y=(finalfactor * cn1peakheight100[ppos])+1.7,
             labels=signif(cn1peakpos100[ppos],3),cex=0.7,srt=0,col="#1C3660AA");
        }
        ppos = c();
        if($only_label_highest_peak){
            ppos = which((cn3peakheight100 == max(cn3peakheight100)) & (cn3peakheight100 > $minimum_labelled_peak_height));
        } else {
            ppos = which(cn3peakheight100 > $minimum_labelled_peak_height);
        }
        if(!(length(ppos) == 0)){
            text(x=cn3peakpos100[ppos],
                 y=(finalfactor * cn3peakheight100[ppos])+1.7,
             labels=signif(cn3peakpos100[ppos],3),cex=0.7,srt=0,col="#F49819AA");
        }
        ppos = c();
        if($only_label_highest_peak){
            ppos = which((cn4peakheight100 == max(cn4peakheight100)) & (cn4peakheight100 > $minimum_labelled_peak_height));
        } else {
            ppos = which(cn4peakheight100 > $minimum_labelled_peak_height);
        }
        if(!(length(ppos) == 0)){
            text(x=cn4peakpos100[ppos],
                 y=(finalfactor * cn4peakheight100[ppos])+1.7,
             labels=signif(cn4peakpos100[ppos],3),cex=0.7,srt=0,col="#E52420AA");
        }
    }

    ppos = c();
    if($only_label_highest_peak){
        ppos = which((cn2peakheight100 == max(cn2peakheight100)) & (cn2peakheight100 > $minimum_labelled_peak_height));
    } else {
        ppos = which(cn2peakheight100 > $minimum_labelled_peak_height);
    }
    if(!(length(ppos) == 0)){
        text(x=cn2peakpos100[ppos],
             y=(finalfactor * cn2peakheight100[ppos])+1.7,
             labels=signif(cn2peakpos100[ppos],3),cex=0.7,srt=0,col="#67B32EAA");
    }

    axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
    mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
    mtext("Kernel Density\n",side=2,cex=0.5,padj=0);
    mtext(genome,adj=0,padj=-5,cex=0.65,side=3);



_END_OF_R_
#-------------------------------------------------
    print R_COMMANDS "$R_command\n";


    #if cn is being plotted
    if((defined($cnvhmm_file) || (defined($cbs_file))) && !($plot_only_CN2)){
        print R_COMMANDS 'drawPlot(z1, cn1minus100x, cn1xchr100, additional_plot_points_cn1, cncircle=1)' . "\n";
        print R_COMMANDS 'drawPlot(z1, cn2100x, cn2xchr100, additional_plot_points_cn2, cncircle=2, num_clusters=num_clusters, output_filename=clustered_data_output_file)' . "\n";
        print R_COMMANDS 'drawPlot(z1, cn3100x, cn3xchr100, additional_plot_points_cn3, cncircle=3)' . "\n";
        print R_COMMANDS 'drawPlot(z1, cn4plus100x, cn4xchr100, additional_plot_points_cn4, cncircle=4)' . "\n";
    } else {
        print R_COMMANDS 'drawPlot(z1, cn2100x, cn2xchr100, additional_plot_points_cn2, num_clusters=num_clusters, output_filename=clustered_data_output_file)' . "\n";
    }
#-------------------------------------------------
    $R_command = <<"_END_OF_R_";
    axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=-1.2);
mtext("Tumor Variant Allele Frequency",adj=0.5,padj=3.2,cex=0.5,side=1);

_END_OF_R_
#-------------------------------------------------
    print R_COMMANDS "$R_command\n";

} elsif ($analysis_type eq 'wgs') {

#-------------------------------------------------
    $R_command = <<"_END_OF_R_";

    #all coverage points plotted
    finalfactor = 25 / maxden;

    plot.default(x=c(1:10),y=c(1:10),ylim=c(0,28),xlim=c(0,100),axes=FALSE, ann=FALSE,col="#00000000",xaxs="i",yaxs="i");
    rect(0, 0, 100, 28, col = "#00000011",border=NA); #plot bg color
    #lines(c(10,100),c(25,25),lty=2,col="black");
    axis(side=2,at=c(0,25),labels=c(0,sprintf("%.3f", maxden)),las=1,cex.axis=0.6,hadj=0.6,lwd=0.5,lwd.ticks=0.5,tck=-0.01);
    lines(den2\$x,(finalfactor * den2factor),col="#67B32EAA",lwd=2);
    if(!($plot_only_CN2)){
        lines(den1\$x,(finalfactor * den1factor),col="#1C3660AA",lwd=2);
        lines(den3\$x,(finalfactor * den3factor),col="#F49819AA",lwd=2);
        lines(den4\$x,(finalfactor * den4factor),col="#E52420AA",lwd=2);
    }

    #legend
    if ($plot_only_CN2) {
        legend(x="topright",lwd=2,legend=c("2 Copies"), col=c("#67B32EAA"),bty="n",cex=0.6, y.intersp=1.25);
    } else { 
        legend(x="topright",lwd=2,legend=c("1 Copy","2 Copies","3 Copies","4 Copies"), col=c("#1C3660AA","#67B32EAA","#F49819AA","#E52420AA"),bty="n",cex=0.6, y.intersp=1.25);
    }

    ppos = c();
    if($only_label_highest_peak){
        ppos = which((cn1peakheight == max(cn1peakheight)) & (cn1peakheight > $minimum_labelled_peak_height));
    } else {
        ppos = which(cn1peakheight > $minimum_labelled_peak_height);
    }
    if(!($plot_only_CN2)){
        if(!(length(ppos) == 0)){    
            text(x=cn1peakpos[ppos], y=(finalfactor * cn1peakheight[ppos])+1.7,
             labels=signif(cn1peakpos[ppos],3),
             cex=0.7, srt=0, col="#1C3660AA");
        }
        ppos = c();
        if($only_label_highest_peak){
            ppos = which((cn3peakheight == max(cn3peakheight)) & (cn3peakheight > $minimum_labelled_peak_height));
        } else {
            ppos = which(cn3peakheight > $minimum_labelled_peak_height);
        }
        if(!(length(ppos) == 0)){    
            text(x=cn3peakpos[ppos],
                 y=(finalfactor * cn3peakheight[ppos])+1.7,
             labels=signif(cn3peakpos[ppos],3),
             cex=0.7,srt=0,col="#F49819AA");
        }
        ppos = c();
        if($only_label_highest_peak){
            ppos = which((cn4peakheight == max(cn4peakheight)) & (cn4peakheight > $minimum_labelled_peak_height));
        } else {
            ppos = which(cn4peakheight > $minimum_labelled_peak_height);
        }
        if(!(length(ppos) == 0)){    
            text(x=cn4peakpos[ppos],
                 y=(finalfactor * cn4peakheight[ppos])+1.7,
             labels=signif(cn4peakpos[ppos],3),
             cex=0.7,srt=0,col="#E52420AA");
        }
    }
    ##print(cn2peakpos);
    ppos = c();
    if($only_label_highest_peak){
        ppos = which((cn2peakheight == max(cn2peakheight)) & (cn2peakheight > $minimum_labelled_peak_height));
        ##cat(max(cn2peakheight),"--",cn2peakheight[ppos],"--",cn2peakpos[ppos])
    } else {
        ppos = which(cn2peakheight > $minimum_labelled_peak_height);
    }
    if(!(length(ppos) == 0)){   
        text(x=cn2peakpos[ppos],
             y=(finalfactor * cn2peakheight[ppos])+1.7,
             labels=signif(cn2peakpos[ppos],3),
             cex=0.7,srt=0,col="#67B32EAA");
    }


    axis(side=3,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=0.6,lwd=0.5,lwd.ticks=0.5,padj=1.4);
    mtext("Tumor Variant Allele Frequency",adj=0.5,padj=-3.1,cex=0.5,side=3);
    mtext("Kernel Density\n",side=2,cex=0.5,padj=0);
    mtext(genome,adj=0,padj=-5,cex=0.65,side=3);


_END_OF_R_
#-------------------------------------------------
    print R_COMMANDS "$R_command\n";

    #if cn is being plotted
    if((defined($cnvhmm_file) || defined($cbs_file)) && !($plot_only_CN2)){
        print R_COMMANDS 'drawPlot(z1, cn1minus, cn1xchr, additional_plot_points_cn1, cncircle=1)' . "\n";
        print R_COMMANDS 'drawPlot(z1, cn2, cn2xchr, additional_plot_points_cn2, cncircle=2, num_clusters=num_clusters, output_filename=clustered_data_output_file)' . "\n";        
        print R_COMMANDS 'drawPlot(z1, cn3, cn3xchr, additional_plot_points_cn3, cncircle=3)' . "\n";
        print R_COMMANDS 'drawPlot(z1, cn4plus, cn4xchr, additional_plot_points_cn4, cncircle=4)' . "\n";
    } else {
        print R_COMMANDS 'drawPlot(z1, cn2, cn2xchr, additional_plot_points_cn2, num_clusters=num_clusters, output_filename=clustered_data_output_file)' . "\n";
    }
}
#-------------------------------------------------
$R_command = <<"_END_OF_R_";
devoff <- dev.off();
q();
_END_OF_R_
#-------------------------------------------------

print R_COMMANDS "$R_command\n";

close R_COMMANDS;

my $cmd = "R --vanilla --slave \< $r_script_output_file";
my $return = Genome::Sys->shellcmd(
    cmd => "$cmd",
    output_files => [$output_image],
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
        my ($chr,$start,$end,$size,$nmarkers,$cn,$adjusted_cn,$cn_normal,$adjusted_cn_normal,$score)=split(/\t/, $line);
        my $pos=$start."_".$end;
        $info_hash{$chr}{$pos}=$adjusted_cn;
    }
    return \%info_hash;
}


sub build_hash_cbs
{
    my ($file)=@_;
    my %info_hash;
    my $fh=new FileHandle($file);
    while(my $line = <$fh>)
    {
        chomp($line);
        unless ($line =~ /^#/){ next;}
        my ($chr,$start,$end,$nmarkers,$adjusted_cn);
        #convert from log2 to abs copy number here 
        $adjusted_cn = (2^$adjusted_cn)*2;
        #round to nearest integer
        $adjusted_cn = sprintf("%.0f", $adjusted_cn);
        my $pos=$start."_".$end;
        $info_hash{$chr}{$pos}=$adjusted_cn;
    }
    return \%info_hash;
}




