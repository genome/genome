package Genome::VariantReporting::Process::CreateReport;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate_pos);

class Genome::VariantReporting::Process::CreateReport {
    is => ['Genome::Process'],

    has_input => [
        input_vcf => {
            is => 'Path',
        },
        variant_type => {
            is => 'Text',
        },
    ],
    has_transient => [
        _plan_file => {
            is => 'Path',
        },
        _translations_file => {
            is => 'Path',
        },
    ],
};

sub result_class {
    return "Genome::VariantReporting::Framework::ReportResult";
}

sub output_directory {
    my $self = shift;
    my $generate_reports_result = $self->results(subclass_name => $self->result_class);
    return $generate_reports_result->output_dir;
}

sub save_plan_file {
    my ($self, $plan_file) = validate_pos(@_, 1, 1);
    return $self->_plan_file($plan_file);
}

sub plan_file {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'plan.yaml');
}

sub save_translations_file {
    my ($self, $translations_file) = validate_pos(@_, 1, 1);
    return $self->_translations_file($translations_file);
}

sub translations_file {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'translations.yaml');
}

sub create_disk_allocation {
    my $self = shift;

    my $rv = $self->SUPER::create_disk_allocation();
    Genome::Sys->copy_file($self->_plan_file, $self->plan_file);
    Genome::Sys->copy_file($self->_translations_file, $self->translations_file);
    $self->debug_message("Saved Plan and Translations file to: %s",
        $rv->absolute_path);
    return $rv;
}



1;
