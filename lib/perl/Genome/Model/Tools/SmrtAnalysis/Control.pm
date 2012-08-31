package Genome::Model::Tools::SmrtAnalysis::Control;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::Control {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        filtered_fofn => {
            is => 'Text',
            doc => 'Specify a regions table for filtering reads.',
        },
        control_reference_directory => {
            is => 'Text',
            doc => 'The directory where the PacBio control reference lives.',
        },
    ],
    has_optional_input => [
        algorithm => {
            is => 'Text',
            doc => 'The algorithm to use for alignment.',
            default_value => 'blasr',
        },
        algorithm_params => {
            is => 'Text',
            doc => 'Pass through options for algorithm (e.g. exonerate)',
            default_value => '-bestn 1',
        },
        nproc => {
            is => 'Number',
            doc => 'Number of processors to use for alignment computation',
            default_value => 4,
        },
        min_accuracy => {
            is => 'Number',
            doc => 'Min accuracy to output a hit',
            default_value => 0.75,
        },
        min_length => {
            is => 'Number',
            doc => 'Min length to output a hit',
            default_value => 50,
        },
        min_z => {
            is => 'Number',
            doc => 'Min Z score to output a hit',
            default_value => 3,
        },
        noise_data => {
            is => 'Text',
            doc => 'noise triplet, .xy or compare XML file containing noise data for estimating a z-score',
            default_value => '-77.27,0.08654,0.00121',
        },
        xml => {
            is => 'Boolean',
            doc => 'Generate XML output.',
            default_value => 0,
        },
        hdf5_mode => {
            is => 'Text',
            doc => "'w' for creating a new hdf5 file, 'a' for appending",
            valid_values => ['w','a'],
            default_value => 'w',
        },
        split_subreads => {
            is => 'Boolean',
            doc => 'Split reads into subreads if subread regions are available.',
            default_value => 0,
        },
        filter_adapter_only => {
            is => 'Boolean',
            doc => 'If specified, do not report adapter-only hits using annotations associated with the reference entry.',
            default_value => 1,
        },
        load_pulses_metrics => {
            default_value => 'QualityValue',
        }
    ],
    has_optional_output => {
        cmp_hdf5_file => { },
        post_control_fofn => { },
    },
    has_optional_param => [
        lsf_queue => {
            default_value => 'workflow',
        },
        lsf_resource => {
            default_value => '',
        },
    ],
};

sub help_brief {
    'This module filters reads based on the minimum readlength and read quality you specify.'
}


sub help_detail {
    return <<EOS 

EOS
}

sub execute {
    my $self = shift;

    my $job_directory = $self->job_directory;
    my $output_directory = $job_directory .'/data';
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }
    my $input_fofn_fh = Genome::Sys->open_file_for_reading($self->input_fofn);
    my $filtered_fofn_fh = Genome::Sys->open_file_for_reading($self->filtered_fofn);
    my @input_fofns;
    my @filtered_fofns;
    while (my $input_line = $input_fofn_fh->getline) {
        chomp($input_line);

        my $filtered_line = $filtered_fofn_fh->getline;
        chomp($filtered_line);

        my (undef, $base_tmp_path) = tempfile('XXXXXX', DIR => $output_directory);
        unlink ($base_tmp_path) || die('Failed to remove temp base path '. $base_tmp_path);
        my $tmp_input_fofn = $base_tmp_path .'_input.fofn';
        my $tmp_input_fh = Genome::Sys->open_file_for_writing($tmp_input_fofn);
        print $tmp_input_fh $input_line ."\n";
        push @input_fofns, $tmp_input_fofn;
        $tmp_input_fh->close;

        my $tmp_filtered_fofn = $base_tmp_path .'_filtered.fofn';
        my $tmp_filtered_fh = Genome::Sys->open_file_for_writing($tmp_filtered_fofn);
        print $tmp_filtered_fh $filtered_line ."\n";
        push @filtered_fofns, $tmp_filtered_fofn;
        $tmp_filtered_fh->close;
    }
    my $movie_summary = $output_directory .'/control_results_by_movie.csv';
    my $control_fofn = $output_directory .'/post_control_regions.fofn';
    my $cmp_hdf5_file = $output_directory .'/control_reads.cmp.h5';
    my %params = (
        input_fofns => \@input_fofns,
        lookup_region_table => '_filtered.fofn',
        data_directory => $output_directory,
        algorithm => $self->algorithm,
        algorithm_params => $self->algorithm_params,
        control_reference_directory => $self->control_reference_directory,
        nproc => $self->nproc,
        min_accuracy => $self->min_accuracy,
        min_length => $self->min_length,
        min_z => $self->min_z,
        noise_data => $self->noise_data,
        xml => $self->xml,
        hdf5_mode => $self->hdf5_mode,
        split_subreads => $self->split_subreads,
        filter_adapter_only => $self->filter_adapter_only,
        load_pulses_metrics => $self->load_pulses_metrics,
        control_fofn => $control_fofn,
        movie_summary => $movie_summary,
        cmp_hdf5_file => $cmp_hdf5_file,
    );
    my $module_path = $self->get_class_object->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;
    my $workflow = Workflow::Operation->create_from_xml($xml_path);
    my @errors = $workflow->validate;
    unless ($workflow->is_valid) {
        die('Errors encountered while validating workflow '. $xml_path ."\n". join("\n", @errors));
    }
    my $output = Workflow::Simple::run_workflow_lsf($xml_path,%params);
    unless (defined $output) {
        my @errors = @Workflow::Simple::ERROR;
        for (@errors) {
            print STDERR $_->error ."\n";
        }
        return;
    }
    $self->cmp_hdf5_file($output->{cmp_hdf5_file});
    $self->post_control_fofn($output->{control_fofn});
    my @files_to_unlink;
    push @files_to_unlink, @filtered_fofns;
    push @files_to_unlink, @input_fofns;
    for my $file (@files_to_unlink) {
        unless (unlink($file)) {
            die('Failed to remove intermediate file '. $file);
        }
    }
    return 1;
}
