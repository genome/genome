package Genome::Model::Tools::SmrtAnalysis::Filter;

use strict;
use warnings;

use Genome;

use File::Temp qw/tempfile/;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::Filter {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        min_length => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum Readlength',
            default_value => 50,
        },
        read_score => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum Read Quality',
            default_value => 0.75,
        },
        read_white_list => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    has_optional_output => [
        filtered_summary_file => { },
        filtered_fofn_file => {},
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
    my @input_fofns;
    while (my $line = $input_fofn_fh->getline) {
        chomp($line);
        my ($tmp_input_fh, $tmp_input_fofn) = tempfile('input-XXXXXX', DIR => $output_directory, SUFFIX => '.fofn');
        print $tmp_input_fh $line ."\n";
        push @input_fofns, $tmp_input_fofn;
        $tmp_input_fh->close;
    }
    
    my $filtered_summary_file = $output_directory .'/filtered_summary.csv';
    if (-e $filtered_summary_file) {
        # Should we just die here instead...
        warn('Removing existing filtered summary file '. $filtered_summary_file);
        unless (unlink($filtered_summary_file)) {
            die('Failed to remove existing filtered summary file '. $filtered_summary_file);
        }
    }
    my $filtered_fofn_file = $output_directory .'/filtered_regions.fofn';
    if (-e $filtered_fofn_file) {
        # Should we just die here instead...
        warn('Removing existing filtered fofn file '. $filtered_fofn_file);
        unless (unlink($filtered_fofn_file)) {
            die('Failed to remove existing filtered fofn file '. $filtered_fofn_file);
        }
    }
    my $filtered_fasta_file = $output_directory .'/filtered_subreads.fa';
    if (-e $filtered_fasta_file) {
        # Should we just die here instead...
        warn('Removing existing filtered fofn file '. $filtered_fasta_file);
        unless (unlink($filtered_fasta_file)) {
            die('Failed to remove existing filtered fofn file '. $filtered_fasta_file);
        }
    }
    my %params = (
        input_fofns => \@input_fofns,
        min_read_length => $self->min_length,
        min_read_score => $self->read_score,
        data_directory => $output_directory,
        filtered_summary_file => $filtered_summary_file,
        filtered_fofn_file => $filtered_fofn_file,
        filtered_fasta_file => $filtered_fasta_file,
        read_white_list => $self->read_white_list,
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
    my @files_to_unlink;
    push @files_to_unlink, @input_fofns;
    for my $file (@files_to_unlink) {
        unless (unlink($file)) {
            die('Failed to remove intermediate file '. $file);
        }
    }
    $self->filtered_summary_file($output->{filtered_summary_file});
    $self->filtered_fofn_file($output->{filtered_fofn_file});
    return 1;
}
