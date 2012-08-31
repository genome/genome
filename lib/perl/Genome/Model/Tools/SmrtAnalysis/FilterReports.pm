package Genome::Model::Tools::SmrtAnalysis::FilterReports;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::FilterReports {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        filtered_summary_file => {
            is => 'Text',
            doc => 'The filtered summary csv file.',
        }
    ],
    has_optional_output => [
        filter_report_xml_file => { },
        adapter_report_xml_file => { },
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
    my $output_directory = $job_directory .'/results';
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }
    my %params = (
        input_fofn => $self->input_fofn,
        results_directory => $output_directory,
        filtered_summary_file => $self->filtered_summary_file,
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
    $self->filter_report_xml_file($output->{filter_report_xml_file});
    $self->adapter_report_xml_file($output->{adapter_report_xml_file});
    return 1;
}
