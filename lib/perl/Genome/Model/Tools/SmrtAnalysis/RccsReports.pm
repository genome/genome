package Genome::Model::Tools::SmrtAnalysis::RccsReports;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::RccsReports {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        rccs_per_base_info_file => {
            is => 'Text',
            doc => 'The per base info file from RCCS.',
        },
    ],
    has_optional_output => [
        rccs_report_xml_file => {},
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
        results_directory => $output_directory,
        rccs_per_base_info_file => $self->rccs_per_base_info_file,
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
    $self->rccs_report_xml_file($output->{rccs_report_xml_file});
    return 1;
}


1;
