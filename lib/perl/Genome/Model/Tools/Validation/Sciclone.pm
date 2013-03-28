package Genome::Model::Tools::Validation::Sciclone;

use strict;
use warnings;
use FileHandle;
use Genome;
use File::Basename;
use FileHandle;

class Genome::Model::Tools::Validation::Sciclone {
    is => 'Command',

    has => [

        variant_files => {
            is => 'Text',
            doc => "comma separated list - files of validated variants with readcounts. 7-column Bam-readcount format - columns: Chr, Start, Ref, Var, RefReads, VarReads, VAF",
            is_optional => 0,
            is_input => 1 ,
        },

        copy_number_files => {
            is => 'Text',
            doc => 'comma separated list - files of segmented copy number calls. If not specified, assumes all variants are CN 2. Expects 5-col format - Chr, St, Sp, NumProbes, SegMean. If you have CNVHMM calls, use "gmt copy-number convert-cnvhmm-output-to-sane-format" to fix them',
            is_optional => 1,
            is_input => 1
        },

        sample_names => {
            is => 'Text',
            doc => "comma separated list - Sample name to be put on graphs",
            is_optional => 1,
            is_input => 1,
        },

        highlight_sex_chrs => {
            is => 'Boolean',
            doc => "Highlight the sex chromosomes (X|Y) on the plot",
            is_optional => 1,
            is_input => 1,
            default => 0
        },

        regions_to_exclude => {
            is => 'Text',
            doc => "comma separated list - regions to exclude (first 3 cols are chr,st,sp). Commonly used for LOH calls",
            is_optional => 1,
            is_input => 1,
            default => 1,
        },

        # positions_to_highlight => {
        #     is => 'Text',
        #     doc => "A tab-delimited file listing variants highlight on the plot. First two columns must be Chr, St",
        #     is_optional => 1,
        #     is_input => 1
        # },


        output_prefix => {
            is => 'Text',
            doc => "prefix for output files",
            is_optional => 0,
            is_input => 1,
            is_output => 1
        },

        skip_if_output_is_present => {
            is => 'Text',
            doc => "Skip if Output is Present",
            is_optional => 1,
            is_input => 1,
            default => 0,
        },

        # minimum_labelled_peak_height => {
        #     is => 'Text',
        #     doc => "only peaks that exceed this height get labelled",
        #     is_optional => 1,
        #     is_input => 1,
        #     default => 0.001
        # },

        # only_label_highest_peak => {
        #     is => 'Boolean',
        #     doc => "only label the highest peak",
        #     is_optional => 1,
        #     default => 0
        # },

        # plot_only_cn2 => {
        #     is => 'Boolean',
        #     doc => "only plot the CN2 data",
        #     is_optional => 1,
        #     default => 0
        # },

        overlay_clusters => {
            is => 'Boolean',
            doc => "overlay information about how the points clustered onto the scatterplot",
            is_optional => 1,
            default => 1,
        },

        tumor_purities => {
            is => 'Integer',
            doc => "comma separated list of tumor purities (between 0 to 100). Will be estimated by the tool if not provided",
            is_optional => 1,
        },

        # component_distribution => {
        #     is => 'String',
        #     doc => "If clustering data, this will be the distribution used for fitting components. Please choose between 'Binomial' or 'Normal'.",
        #     is_optional => 1,
        #     default => "Binomial"
        # },

        minimum_depth => {
            is => 'Integer',
            doc => "Plot/Cluster only using variants that have at least this many reads. 100 is a reasonable default for capture data. If only wgs data is available, you'll need to lower this value.",
            is_optional => 1,
            default => 100,
        },

        cn_calls_are_log2 => {
            is => 'Boolean',
            doc => "copy number calls are in log 2 format (default is absolute CN - 1, 2, 3, etc)",
            is_optional => 1,
            default => 0,
        },

        do_clustering => {
            is => 'Boolean',
            doc => "if true, clusters the data. if false, just creates a plot (saving time)",
            is_optional => 1,
            default => 1,
        },

        # label_highlighted_points => {
        #     is => 'Boolean',
        #     doc => "if true, assumes that the positions-to-highlight file has a third column, containing names for the points to be highlighted. Sets plotOnlyCN2 to true and adds numbered labels and a legend for the highlighted points",
        #     is_optional => 1,
        #     default => 0,
        # },


        ],
};


sub help_brief {
    "Plot CN-separated SNV density plots with optional clustering"
}

sub help_synopsis {
    return <<EOS
        Inputs of variant readcounts and copy-number segmentation data, Output of pdf plots.
EXAMPLE:	gmt validation clonality-plot --variant-file snvs.txt,snvs2.txt --output-prefix clonality --sample-names 'Sample1,Sample2' --copy_number_files segs.paired.dat,segs2.paired.dat
EOS
}

sub help_detail {
    return <<EOS
This tool can be used to plot the CN-separated SNV density plots that are known at TGI as 'clonality plots'. Can be used for WGS or Capture data, but best results will be had with greater than 100x read depth.
EOS
}

sub execute {
    my $self = shift;

    ##inputs##
    my $variant_files = $self->variant_files;
    my $cn_files = $self->copy_number_files;
    my $sample_names = $self->sample_names;
    my $highlight_sex_chrs = $self->highlight_sex_chrs;
    #my $positions_to_highlight_file = $self->positions_to_highlight;
    my $tumor_purities = $self->tumor_purities;
    #my $component_distribution = $self->component_distribution;
    #my $maximum_clusters_to_test = $self->max_clusters_to_test;
    my $cn_calls_are_log2 = $self->cn_calls_are_log2;
    #my $use_sex_chrs = $self->use_sex_chrs;
    my $regions_to_exclude = $self->regions_to_exclude;

    ##outputs##
    my $output_prefix = $self->output_prefix;
    my $r_script_output_file = "$output_prefix.R";

    ##options##
    my $skip_if_output_is_present = $self->skip_if_output_is_present;
    #my $minimum_labelled_peak_height = $self->minimum_labelled_peak_height;
    #my $only_label_highest_peak = $self->only_label_highest_peak;
    #my $plot_only_cn2 = $self->plot_only_cn2;
    my $overlay_clusters = $self->overlay_clusters;
    my $minimum_depth = $self->minimum_depth;
    my $do_clustering = $self->do_clustering;
    #my $label_highlighted_points = $self->label_highlighted_points;


    # if(($label_highlighted_points) && !(defined($positions_to_highlight_file))){
    #     die("ERROR: if label-highlighted-points is true, a positions-to-highlight file must be provided");
    # }

    my $rfile;
    open($rfile, ">$r_script_output_file") || die "Can't open R file for writing.\n";

    # write out the r commands
    print $rfile "library(sciClone)\n";

    ##TODO - handle headers
    
    #read in the variant files
    my @variantFiles = split(",",$variant_files);
    my @variantVars;
    my $i=0;
    for($i=0;$i<@variantFiles;$i++){        
        my $var = "v$i";
        push(@variantVars,$var);
        # read in the file (and convert varscan, if necessary)
        print $rfile "$var = " . 'read.table("' . $variantFiles[$i] .  '")' . "\n";
        print $rfile "$var = $var" . '[,c(1,2,5,6,7)]' . "\n";
    }

    #read in the cn files
    my @cnVars;
    if(defined($cn_files)){
        my @cnFiles = split(",",$cn_files);
        my $i=0;
        for($i=0;$i<@cnFiles;$i++){        
            my $var = "cn$i";
            push(@cnVars,$var);
            # read in the file (and convert varscan, if necessary)
            print $rfile "$var = " . 'read.table("' . $cnFiles[$i] .  '")' . "\n";
            print $rfile "$var = $var" . '[,c(1,2,3,5)]' . "\n";
        }
    }

    my @regVars;
    if(defined($regions_to_exclude)){
        my @regFiles =  split(",",$regions_to_exclude);
        my $i=0;
        for($i=0;$i<@regFiles;$i++){        
            my $var = "reg$i";
            push(@regVars,$var);
            # read in the file (and convert varscan, if necessary)
            print $rfile "$var = " . 'read.table("' . $regFiles[$i] .  '")' . "\n";
            print $rfile "$var = $var" . '[,c(1,2,3)]' . "\n";
        }        
    }

    my @sampleNames = split(",",$sample_names);
    my $sampleNames = '"' . join('","',@sampleNames) . '"';

    #print out the clonecaller command, one piece at a time
    my $cmd = "sciClone(outputPrefix=$output_prefix";
    $cmd = $cmd . ", vafs=list(" . join(",",@variantVars) . ")";
    $cmd = $cmd . ", sampleNames=c(" . $sampleNames . ")";

    if(defined($cn_files)){
        $cmd = $cmd . ", copyNumberCalls=list(" . join(",",@cnVars) . ")";
    }

    if(defined($regions_to_exclude)){
        $cmd = $cmd . ", regionsToExclude=list(" . join(",",@regVars) . ")";
    }

    $cmd = $cmd . ", minimumDepth=$minimum_depth";

    if(defined($tumor_purities)){
        print "tp: " . $tumor_purities . "\n";
        $cmd = $cmd . ", purity=c($tumor_purities)";
    }

    if($overlay_clusters){
        $cmd = $cmd . ", overlayClusters=TRUE";
    }

    unless($highlight_sex_chrs){
        $cmd = $cmd . ", highlightSexChrs=FALSE";
    }

    if($cn_calls_are_log2){
        $cmd = $cmd . ", cnCallsAreLog2=TRUE";
    }

    unless($do_clustering){
        $cmd = $cmd . ", doClustering=FALSE";
    }

    print $rfile $cmd . ")\n";
    close $rfile;

    #now actually run the R script
    # my $rcmd = "R --vanilla --slave \< $r_script_output_file";
    my $rcmd = "Rscript $r_script_output_file";
    my $return_value = Genome::Sys->shellcmd(
        cmd => "$rcmd",
        output_files => ["$output_prefix.clusters"],
        skip_if_output_is_present => $skip_if_output_is_present,
        );
    unless($return_value) {
        $self->error_message("Failed to execute: Returned $return_value");
        die $self->error_message;
    }
    return $return_value;
}
