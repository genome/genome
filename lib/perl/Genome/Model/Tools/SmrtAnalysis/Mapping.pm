package Genome::Model::Tools::SmrtAnalysis::Mapping;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::Mapping {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        post_control_fofn => {
            is => 'Text',
            doc => 'Specify a regions table for filtering reads.',
        },
        reference_directory => {
            is => 'Text',
            doc => 'The directory where the PacBio reference lives.',
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
            default_value => '-bestn 1 -minMatch 8 -minPctIdentity 70.0',
        },
        nproc => {
            is => 'Number',
            doc => 'Number of processors to use for alignment computation',
            default_value => 4,
        },
        use_ccs => {
            is => 'Text',
            valid_values => ['fullpass','allpass','denovo'],
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
        pulse_metrics => {
            default_value => 'QualityValue,InsertionQV,DeletionQV,IPD,PulseWidth',
        }
    ],
    has_optional_output => [
        cmp_hdf5_file => { },
        alignment_summary_gff => { },
    ],
    has_optional_param => [
        lsf_queue => {
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW},
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
    my $post_control_fofn_fh = Genome::Sys->open_file_for_reading($self->post_control_fofn);
    my @input_fofns;
    my @post_control_fofns;
    while (my $input_line = $input_fofn_fh->getline) {
        chomp($input_line);

        my $post_control_line = $post_control_fofn_fh->getline;
        chomp($post_control_line);

        my (undef, $base_tmp_path) = tempfile('XXXXXX', DIR => $output_directory);
        unlink ($base_tmp_path) || die('Failed to remove temp base path '. $base_tmp_path);
        my $tmp_input_fofn = $base_tmp_path .'_input.fofn';
        my $tmp_input_fh = Genome::Sys->open_file_for_writing($tmp_input_fofn);
        print $tmp_input_fh $input_line ."\n";
        push @input_fofns, $tmp_input_fofn;
        $tmp_input_fh->close;

        my $tmp_post_control_fofn = $base_tmp_path .'_post_control.fofn';
        my $tmp_post_control_fh = Genome::Sys->open_file_for_writing($tmp_post_control_fofn);
        print $tmp_post_control_fh $post_control_line ."\n";
        push @post_control_fofns, $tmp_post_control_fofn;
        $tmp_post_control_fh->close;
    }
    my @suffix_arrays = glob($self->reference_directory .'/sequence/*.fasta.sa');
    my $algorithm_params = $self->algorithm_params;
    if (@suffix_arrays) {
        unless (scalar(@suffix_arrays) == 1) {
            warn('Please FIX!  More than one suffix array found.');
            die(Data::Dumper::Dumper(@suffix_arrays));
        }
        $algorithm_params .= ' -sa '. $suffix_arrays[0];
    }
    $self->algorithm_params($algorithm_params);

    my $cmp_hdf5_file = $output_directory .'/aligned_reads.cmp.h5';
    my $alignment_summary_gff = $output_directory .'/alignment_summary.gff';
    my $coverage_bed_file = $output_directory .'/coverage.bed';
    my %params = (
        input_fofns => \@input_fofns,
        lookup_region_table => '_post_control.fofn',
        data_directory => $output_directory,
        algorithm => $self->algorithm,
        algorithm_params => $self->algorithm_params,
        reference_directory => $self->reference_directory,
        nproc => $self->nproc,
        min_accuracy => $self->min_accuracy,
        min_length => $self->min_length,
        xml => $self->xml,
        hdf5_mode => $self->hdf5_mode,
        cmp_hdf5_file => $cmp_hdf5_file,
        alignment_summary_gff => $alignment_summary_gff,
        coverage_bed_file => $coverage_bed_file,
        gff_to_bed_purpose => 'coverage',
        gff_to_bed_name => 'meanCoverage',
        gff_to_bed_description => 'Mean coverage of genome in fixed interval regions',
        pulse_metrics => $self->pulse_metrics,
    );
    if ($self->use_ccs) {
        $params{'use_ccs'} = $self->use_ccs;
    }
    my $module_path = $self->__meta__->module_path;
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
    $self->alignment_summary_gff($output->{alignment_summary_gff});
    my @files_to_unlink;
    push @files_to_unlink, @post_control_fofns;
    push @files_to_unlink, @input_fofns;
    for my $file (@files_to_unlink) {
        unless (unlink($file)) {
            die('Failed to remove intermediate file '. $file);
        }
    }
    return 1;
}
