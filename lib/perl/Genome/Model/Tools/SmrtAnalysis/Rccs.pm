package Genome::Model::Tools::SmrtAnalysis::Rccs;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::Rccs {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            is => 'Text',
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        cmp_hdf5_file => {
            is => 'Text',
            doc => 'An aligned reads cmp.h5 format file.',
        },
        reference_directory => {
            is => 'Text',
            doc => 'The directory where the PacBio reference lives.',
        },
    ],
    has_optional_input => [
        min_subread_accuracy => {
            is => 'Number',
            default_value => 0.75,
        },
        min_fraction => {
            is => 'Number',
            default_value => 0,
        },
        file_prefix => {
            is => 'Text',
            default_value => 'RCCS',
        },
    ],
    has_optional => [
    ],
    has_optional_output => [
        frequency_count_gff_file => { },
        consensus_fastq_file => { },
        consensus_fasta_file => { },
        rccs_per_base_info_file => { },
        rccs_bam_file => { },
    ],
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
    ''
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
    my $raw_read_fasta_file = $output_directory .'/read_sequences.fa';
    my %params = (
        input_fofn => $self->input_fofn,
        raw_read_fasta_file => $raw_read_fasta_file,
        reference_directory => $self->reference_directory,
        data_directory => $output_directory,
        file_prefix => $self->file_prefix,
        cmp_hdf5_file => $self->cmp_hdf5_file,
        min_subread_accuracy => $self->min_subread_accuracy,
        min_fraction => $self->min_fraction,
        split_subreads => 0,
        trim_by_region => 0,
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
    $self->frequency_count_gff_file($output->{frequency_count_gff_file});
    $self->consensus_fastq_file($output->{consensus_fastq_file});
    $self->consensus_fasta_file($output->{consensus_fasta_file});
    $self->rccs_per_base_info_file($output->{rccs_per_base_info_file});
    $self->rccs_bam_file($output->{rccs_bam_file});
    return 1;
}


1;
