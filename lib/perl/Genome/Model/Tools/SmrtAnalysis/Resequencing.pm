package Genome::Model::Tools::SmrtAnalysis::Resequencing;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

class Genome::Model::Tools::SmrtAnalysis::Resequencing {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        input_fofn => {
            doc => 'Eeither pls.h5 or bas.h5 file fofn',
        },
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        control_reference_directory => {
            is => 'Text',
            doc => 'The directory where the PacBio control reference lives.',
        },
        reference_directory => {
            is => 'Text',
            doc => 'The directory where the PacBio reference lives.',
        },
    ],
    has_optional_input => [
        min_length => {
            is => 'Number',
            doc => 'Minimum Readlength',
            default_value => 50,
        },
        read_score => {
            is => 'Number',
            doc => 'Minimum Read Quality',
            default_value => 0.75,
        },
        use_ccs => {
            is => 'Text',
            doc => 'Map the ccsSequence to the genome first, then align subreads to the interval that the CCS read mapped to.  fullpass only aligns subreads that span the length of the template.  Specifying allpass maps all subreads.',
            valid_values => ['fullpass','allpass','denovo'],
        },
    ],
};


sub execute {
    my $self = shift;
    my %params = (
        input_fofn => $self->input_fofn,
        min_length => $self->min_length,
        read_score => $self->read_score,
        job_directory => $self->job_directory,
        reference_directory => $self->reference_directory,
        control_reference_directory => $self->control_reference_directory,
    );
    if ($self->use_ccs) {
        $params{'use_ccs'} = $self->use_ccs;
    }
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
    return 1;
}


1;
