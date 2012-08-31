package Genome::Model::Tools::BioSamtools::CoverageStatsV2;

use strict;
use warnings;

use Genome;
use Workflow::Simple;

my $DEFAULT_MINIMUM_DEPTHS = '1,5,10,15,20';
my $DEFAULT_WINGSPAN_VALUES = '0,200,500';
my $DEFAULT_LOG_DIRECTORY = '/gsc/var/log/genome/coverage_stats';

class Genome::Model::Tools::BioSamtools::CoverageStatsV2 {
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
        genome_file => {
            is => 'Text',
            doc => 'Required if winspan_values is defined greater than zero.  Used to add wingspan via slopBed.',
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
            default_value => $DEFAULT_LOG_DIRECTORY,
            is_optional => 1,
        },
    ],
    has_output => [
        stats_files => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of stats files produced by the workflow'},
        alignment_summaries => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of alignment summaries produced by the workflow'},
        stats_summaries => { is => 'Array', is_many => 1, is_optional => 1, doc => 'a list of stats summaries produced by the workflow'},
    ]
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }

    unless ($] > 5.010) {
        die 'Subcommands run by '. __PACKAGE__ .' require perl 5.10 or greater! Consider using \'/usr/bin/perl `which gmt`\' instead of \'gmt\'.';
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

    my ($bed_basename,$bed_dirname,$bed_suffix) = File::Basename::fileparse($self->bed_file,qw/\.bed/);
    
    my @wingspans = split(',',$self->wingspan_values);
    my @bed_files;
    for my $wingspan (@wingspans) {
        my $wingspan_bed_file = $self->output_directory .'/'. $bed_basename .'-wingspan_'.$wingspan .'.bed';
        if ($wingspan > 0) {
            unless ($self->genome_file) {
                die('In order to add the wingspan via Genome::Model::ToolsBedTool::Slop we need a genome_file defined!');
            }
            my $tmp_slop_bed_file = Genome::Sys->create_temp_file_path();
            my %slop_bed_params = (
                input_file => $self->bed_file,
                genome_file => $self->genome_file,
                output_file => $tmp_slop_bed_file,
                both => $wingspan,
            );
            unless (Genome::Model::Tools::BedTools::Slop->execute(%slop_bed_params)) {
                die('Failed to run SlopBed with params : '. Data::Dumper::Dumper(%slop_bed_params));
            }
            unless (Genome::Model::Tools::BedTools::Merge->execute(
                input_file => $tmp_slop_bed_file,
                output_file => $wingspan_bed_file,
                report_names => 1,
            )) {
                die('Failed to run MergeBed after wingspan slop added!');
            }
        } else {
            File::Copy::copy($self->bed_file,$wingspan_bed_file);
        }
        push @bed_files, $wingspan_bed_file;
    }
    
    my $module_path = $self->__meta__->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;

    my $workflow = Workflow::Operation->create_from_xml($xml_path);

    unless (-d $self->log_directory) {
        unless (Genome::Sys->create_directory($self->log_directory)) {
            die('Failed to create output_directory: '. $self->log_directory);
        }
    }
    $workflow->log_dir($self->log_directory);

    my $output = run_workflow_lsf(
        $workflow,
        bed_files => \@bed_files,
        bam_file => $self->bam_file,
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


sub setup_workflow_operation {
    my $self = shift;
    my %params = @_;

    my $workflow = delete($params{'workflow'});
    my $name = delete($params{'name'});
    my $class_name = delete($params{'class_name'});
    my $input_properties = delete($params{'input_properties'});
    my $output_properties = delete($params{'output_properties'});
    my $parallel_by = delete($params{'parallel_by'});
    
    my $input_connector = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    my $operation = $workflow->add_operation(
        name => $name,
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => $class_name,
        ),
    );
    for my $left_property (keys %{$input_properties}) {
        my $right_property = $input_properties->{$left_property};
        $workflow->add_link(
            left_operation => $input_connector,
            left_property => $left_property,
            right_operation => $operation,
            right_property => $right_property,
        );
    }
    if ($parallel_by) {
        $operation->parallel_by($parallel_by);
    }
    for my $left_property (keys %{$output_properties}) {
        my $right_property = $output_properties->{$left_property};
        $workflow->add_link(
            left_operation => $operation,
            left_property => $left_property,
            right_operation => $output_connector,
            right_property => $right_property,
        );
    }
    return $operation;
}


1;
