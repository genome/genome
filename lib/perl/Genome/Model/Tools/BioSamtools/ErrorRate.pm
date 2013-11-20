package Genome::Model::Tools::BioSamtools::ErrorRate;

use strict;
use warnings;

use Genome;
use Statistics::R;
use Cwd;
use Bio::DB::Sam::Constants;

my $DEFAULT = '0.9';

class Genome::Model::Tools::BioSamtools::ErrorRate {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => {
            is  => 'Text',
            doc => 'A BAM format file of alignment data'
        },
        output_file => {
            is  => 'Text',
            doc => 'A file path to store tab separated value output.  The file extension must be .tsv in pileup mode',
        },
        version => {
            is => 'Version',
            is_optional => 1,
            default_value => $DEFAULT,
            doc => 'version of C pileup tool',
        },
    ],
};


my %versions = (
    0.6 => '/gsc/pkg/bio/bam-errorrate/0.6/bam-errorrate',
    0.7 => '/gsc/pkg/bio/bam-errorrate/0.7/bam-errorrate0.7',
    0.8 => '/gsc/pkg/bio/bam-errorrate/0.8/bam-errorrate0.8',
    0.9 => '/gsc/pkg/bio/bam-errorrate/0.9/bam-errorrate0.9',
    #C util is running out of iferguson home directory until it has been tested and can be deployed to the blades
);


sub default_errorrate_version {
    die "default error_rate version: $DEFAULT is not valid" unless $versions{$DEFAULT};
    return $DEFAULT;
}


sub execute {
    my $self = shift;

    my ($basename, $dirname, $suffix) = File::Basename::fileparse($self->output_file, qw/.tsv/);
    unless ($suffix eq '.tsv') {
        die 'Unable to parse output *.tsv file name: '. $self->output_file;
    }

    my $bam_file    = $self->bam_file;
    my $output_file = $self->output_file;
    my $version     = $self->version;

    my $tool_path   = $versions{$version};
    unless ($tool_path and -s $tool_path and -x $tool_path) {
        die $self->error_message('bam-errorrate version'.$version.' is not valid tool');
    }

    my $cmd = "samtools view $bam_file | $tool_path > $output_file";
    
    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => [ $bam_file ],
        output_files => [ $output_file ],
        skip_if_output_is_present => 0,
    );

    $self->run_r_script;
    return 1;
}


sub run_r_script {
    my $self     = shift;
    my $out_file = $self->output_file;

    my ($basename, $dirname, $suffix) = File::Basename::fileparse($out_file, qw/.tsv/);
    my %read_ends;
    my $confirmed_header = 0;

    my $tsv_fh = Genome::Sys->open_file_for_reading($out_file) 
        or die "failed to open $out_file for reading\n";

    while (my $line = $tsv_fh->getline) {
        next if $line =~ /^#/; # skip comments
        chomp $line;
        my @fields = split /\s+/, $line;

        if (not $confirmed_header) {
            if ($fields[0] ne 'read_end') {
                die sprintf "Unable to parse header in '%s'. Expected first field to be read_end, but it was '%s'.",
                    $out_file,
                    $fields[0];
            }
            $confirmed_header = 1;
            next;
        }

        if ($fields[0] !~ /^[012]$/) {
            die sprintf "Read end was '%s', which was not 0, 1, or 2.", $fields[0];
        }
        $read_ends{$fields[0]} = 1;
    }
    $tsv_fh->close;

    my $r_library = $self->__meta__->module_path;
    $r_library =~ s/\.pm/\.R/;
    $self->status_message('R LIBRARY: '. $r_library);
    my $tempdir = Genome::Sys->create_temp_directory();

    my $cwd = getcwd();

    for my $read_end (keys %read_ends) {
        my $input_file = $out_file;

        my $count_plot_file = $dirname . $basename .'_counts_'. $read_end .'.png';
        my $rate_plot_file  = $dirname . $basename .'_rates_'. $read_end .'.png';
        my $rate_dist_file  = $dirname . $basename .'_rate_distribution_'. $read_end .'.png';

        # The rate plot is all that is necessary at this time
        #my $r_cmd = "generatePlots('$input_file','$read_end','$count_plot_file','$rate_plot_file','$rate_dist_file')";

        my $r_script    = $tempdir .'/'. $basename .'_error_rate_'. $read_end .'.R';
        my $r_script_fh = Genome::Sys->open_file_for_writing($r_script);

        print $r_script_fh "source('$r_library')\n";
        print $r_script_fh "fullFile <- readTable('$input_file')\n";
        print $r_script_fh "readEndErrorRate <- getReadEnd(fullFile,'$read_end')\n";
        print $r_script_fh "positionErrorRate <- getPositionErrorRate(readEndErrorRate)\n";
        print $r_script_fh "makeRatePlot(positionErrorRate,'$rate_plot_file','$read_end')\n";
        #print $r_script_fh "$r_cmd\n";
        $r_script_fh->close;

        my $cmd = 'Rscript '. $r_script;
        Genome::Sys->shellcmd(
            cmd => $cmd,
        );
        unlink $r_script;
    }
    chdir $cwd;
    return 1;
}


1;
