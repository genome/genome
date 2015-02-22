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
        report_names => {
            is => 'Text',
            is_many => 1,
        }
    ],
    has_transient_optional => [
        _plan_file => {
            is => 'Path',
        },
        _translations_file => {
            is => 'Path',
        },
    ],
};

sub symlink_results {
    my $self = shift;
    my $destination = shift;

    Genome::Sys->create_directory($destination);
    for my $report_name ($self->report_names) {
        my $specific_destination = File::Spec->join($destination, $report_name);
        Genome::Sys->create_directory($specific_destination);
        Genome::Sys->symlink_directory($self->output_directory($report_name),
            $specific_destination);
    }

    return 1;
}

sub output_directory {
    my ($self, $report_name) = validate_pos(@_, 1, 1);

    my $result = $self->result_with_label($report_name);
    return $result->output_dir;
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
