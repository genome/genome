package Genome::Model::Tools::BioSamtools::CoverageStats;

use strict;
use warnings;

use Genome;
use Workflow::Simple;

my $DEFAULT_MINIMUM_DEPTHS = '1,5,10,15,20';
my $DEFAULT_WINGSPAN_VALUES = '0,200,500';

class Genome::Model::Tools::BioSamtools::CoverageStats {
    is => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'A path to a BED format file of regions of interest',
        },
        bam_file => {
            is => 'Text',
            doc => 'A path to a BAM format file of aligned capture reads',
        },
        minimum_depths => {
            is => 'Text',
            doc => 'A comma separated list of minimum depths to evaluate coverage',
            default_value => $DEFAULT_MINIMUM_DEPTHS,
            is_optional => 1,
        },
        wingspan_values => {
            is => 'Text',
            doc => 'A comma separated list of wingspan values to add to each region of interest',
            default_value => $DEFAULT_WINGSPAN_VALUES,
            is_optional => 1,
        },
        minimum_base_quality => {
            is => 'Text',
            doc => 'A minimum base quality to consider in coverage assesment',
            default_value => 0,
            is_optional => 1,
        },
        minimum_mapping_quality => {
            is => 'Text',
            doc => 'A minimum mapping quality to consider in coverage assesment',
            default_value => 0,
            is_optional => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'The output directory to generate coverage stats',
        },
        log_directory => {
            is => 'Text',
            doc => 'The directory to store the workflow output and error logs',
            is_optional => 1,
        },
    ],
    has_output => [
        stats_files => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of stats files produced by the workflow'},
        alignment_summaries => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of alignment summaries produced by the workflow'},
        alignment_summaryv2s => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of alignment summaries (new algorithm) produced by the workflow'},
        stats_summaries => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of stats summaries produced by the workflow'},
    ]
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }

    unless ($] > 5.010) {
        die 'Subcommands run by '. __PACKAGE__ .' require perl 5.10 or greater!';
    }
    return $self;
}

sub execute {
    my $self = shift;
    unless (-d $self->output_directory) {
        unless (Genome::Sys->create_directory($self->output_directory)) {
            die('Failed to create output_directory: '. $self->output_directory);
        }
    }
    my $module_path = $self->__meta__->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;
    my @wingspans = split(',',$self->wingspan_values);

    my $workflow = Workflow::Operation->create_from_xml($xml_path);

    if($self->log_directory) {
        unless (-d $self->log_directory) {
            unless (Genome::Sys->create_directory($self->log_directory)) {
                die('Failed to create output_directory: '. $self->log_directory);
            }
        }

        $workflow->log_dir($self->log_directory);
    }

    my $output = run_workflow_lsf($workflow,
                                  bed_file => $self->bed_file,
                                  bam_file => $self->bam_file,
                                  wingspan => \@wingspans,
                                  minimum_depth => $self->minimum_depths,
                                  output_directory => $self->output_directory,
                                  minimum_base_quality => $self->minimum_base_quality,
                                  minimum_mapping_quality => $self->minimum_mapping_quality,
                              );
    unless (defined $output) {
        my @error;
        for (@Workflow::Simple::ERROR) {
            push @error, $_->error;
        }
        $self->error_message(join("\n",@error));
        die($self->error_message);
    }
    my $alignment_summaries = $output->{alignment_summaries};
    unless (scalar(@$alignment_summaries) == scalar(@wingspans)) {
        die('Incorrect number of alignment summaries!');
    }
    $self->alignment_summaries($alignment_summaries);

    my $alignment_summaryv2s = $output->{alignment_summaryv2s};
    unless (scalar(@$alignment_summaryv2s) == scalar(@wingspans)) {
        die('Incorrect number of alignment summaries V2s!');
    }
    $self->alignment_summaryv2s($alignment_summaryv2s);

    my $stats_files_array_ref = $output->{stats_files};
    my @stats_files = @{$stats_files_array_ref};
    unless (scalar(@stats_files) == scalar(@wingspans)) {
        die('Incorrect number of wingspan iterations for stats files!');
    }
    $self->stats_files(\@stats_files);
    my $stats_sum_array_ref = $output->{stats_summaries};
    my @stats_sums = @{$stats_sum_array_ref};
    unless (scalar(@stats_sums) == scalar(@wingspans)) {
        die('Incorrect number of wingspan iterations for stats summariess!');
    }
    $self->stats_summaries(\@stats_sums);
    return 1;
}

1;
