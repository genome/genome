package Genome::Model::Tools::Build::GcBias;

use strict;
use warnings;
use Genome;
use IO::File;
use Statistics::Descriptive;

class Genome::Model::Tools::Build::GcBias {
    is => 'Command',
    has => [
    build_id => {
        type => 'String',
        is_optional => 1,
        doc => "Build-id of the build whose .bam file you would like to calculate a gc-bias graph for.",
    },
    bam_file => {
        type => 'String',
        is_optional => 1,
        doc => "If you do not have a build, input the .bam file, avg read-length, and avg coverage across the genome for this tool. Use full path!",
    },
    mean_read_length => {
        type => 'Number',
        is_optional => 1,
        doc => "Average read length of all lanes included in the .bam file. Input this value when not using a build-id.",
    },
    bam_window_output => {
        type => 'String',
        is_optional => 0,
        doc => "Bam-window output file. Expect a file size of ~46M. Use full path!",
    },
    gc_bias_report_output => {
        type => 'String',
        is_optional => 0,
        doc => "Output file for writing a summary of findings. Use full path!",
    },
    ]
};

sub help_brief {
    "Calculate 10-bin gc bias coverage graph"
}

sub help_detail {
    "RUN ON 64 BIT BLADE. This script takes in a build-id or just a .bam file (plus mean read-length associated with reads in that .bam file) and calculates a coverage value and a count of total reads falling into each of 10 gc-content bins. Bins are divided in terms of gc-content percentage, 0-10, 11-20, and so on. Coverage is reported first non-normalized, and then normalized based on the mean coverage across the entire genome.Read counts are reported in both fashions as well. An offshoot output file for this script is the bam-window output, and also a second file containing the gc-content for each window represented in the bam-window output. This second file is placed next to the bam-window output, having the same filename along with the suffix \".withgc\"."
}

sub execute {

    my $self = shift;

    #test architecture to make sure bam-window program can run (req. 64-bit)
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }

    #declare optional arguments
    my $bam_file;
    my $mean_read_length;
    my $mean_coverage;

    #check for build-id input
    my $build_id = $self->build_id;

    if ($build_id) {

        #get build and model
        my $build = Genome::Model::Build->get($build_id);
        my $model = $build->model;
        unless ($build) {
            $self->error_message("Was unable to find this build: $build_id");
            return;
        }
        unless ($model) {
            $self->error_message("Was unable to find model for this build: $build_id");
            return;
        }

        #get whole rmdup bam file
        $bam_file = $build->whole_rmdup_bam_file;

        #calculate average read length of instrument-data
        my @instrument_data = $model->instrument_data;
        my $read_lengths = Statistics::Descriptive::Full->new();
        for my $lane (@instrument_data) {
            my $rl = $lane->read_length;
            $read_lengths->add_data($rl);
        }
        $mean_read_length = $read_lengths->mean();

        #gather mean coverage across entire genome for build
        my %metrics = map {$_->name => $_->value} $build->metrics;
        $mean_coverage = $metrics{'haploid_coverage'};
    }

    else {

        #use optional arguments for these parameters instead
        $bam_file = $self->bam_file;
        $mean_read_length = $self->mean_read_length;
        unless ($bam_file) {
            $self->error_message("Either build_id or bam_file parameters must be input.");
            return;
        }
        unless ($mean_read_length) {
            $self->error_message("avg_read_length parameter is required when using bam_file parameter.");
            return;
        }
    }

    #run bam-window on bam file
    my $bam_window_output_file = $self->bam_window_output;
    `$ENV{GENOME_SW}/bamwindow/bamwindow-v0.1/bam-window -s -w 1000 -q 1 $bam_file > $bam_window_output_file`;

    #gather gc-content information
    ######## NOTE: this part will FAIL if gc content files are moved ############
    my %gcdata;
    my $gcdir = '/gscmnt/sata407/info/medseq/ndees/GC_and_N_content_files/';
    opendir(GC,$gcdir) or die "Can't open gc-content directory: $gcdir\n";
    while ( defined ( my $file=readdir(GC) ) ) {
        next if ($file=~m/^\./ || $file=~m/whole/);
        $file = $gcdir . $file;
        my $fh = new IO::File $file,"r";
        while (my $line = $fh->getline) {
            my ($chr,$pos,$gc,$n) = split /\t/,$line;
            $gcdata{$chr}{$pos}{gc}=$gc;
        }
    }
    close(GC);

    #fill a hash with info from the bam-window output to be combined with gc-content info
    my %covdata;
    my $bw_fh = new IO::File $bam_window_output_file,"r";
    while (my $line = $bw_fh->getline) {
        next if $line =~ /chr/;
        chomp $line;
        my ($chr,$pos,$cov) = split /\t/,$line;
        $pos = $pos + 999; #doing this for bamwindow's forward-looking positions
        if (exists($gcdata{$chr}{$pos}{gc})) {
            $covdata{$chr}{$pos}{cov} = $cov;
        }
    }
    $bw_fh->close;

    #print coverage plus gc-content info into ".withgc" file
    my $gc_file = $bam_window_output_file . ".withgc";
    my $gc_fh = new IO::File $gc_file,"w";
    print $gc_fh "chr\tpos\tcoverage\tgc_content\n";
    for my $chr (sort keys %covdata) {
        my @positions = keys %{$covdata{$chr}};
        for my $pos (sort { $a <=> $b } @positions) {
            print $gc_fh "$chr\t$pos\t$covdata{$chr}{$pos}{cov}\t$gcdata{$chr}{$pos}{gc}\n";
        }
    }
    $gc_fh->close;

    #call R script to calculate normalized values per gc-content bin
    my $gc_bias_rpt_file = $self->gc_bias_report_output;
    my $rlibrary = "gcbias_lib.R";
    my $cmd;
    #R function: gc_cov_bins = function (in.file,out.file,mean_read_length,mean_coverage=NULL)
    if ($mean_coverage) {
        $cmd = "gc_cov_bins(in.file='$gc_file',out.file='$gc_bias_rpt_file',mean_read_length='$mean_read_length',mean_coverage='$mean_coverage');";
    }
    else {
        $cmd = "gc_cov_bins(in.file='$gc_file',out.file='$gc_bias_rpt_file',mean_read_length='$mean_read_length');";
    }
    my $cmd_Rcall = Genome::Model::Tools::R::CallR->create(command=>$cmd,library=>$rlibrary);
    $cmd_Rcall->execute;

    return 1;
}
1;
