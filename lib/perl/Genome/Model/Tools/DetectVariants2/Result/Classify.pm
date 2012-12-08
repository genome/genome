package Genome::Model::Tools::DetectVariants2::Result::Classify;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Classify {
    is => ['Genome::Model::Tools::DetectVariants2::Result::Base'],
    is_abstract => 1,
    has_input =>[
        prior_result_id => {
            is => 'Text',
            doc => 'ID of the results to be classified',
        },
    ],
    has_param => [
        classifier_version => {
            is => 'Text',
            doc => 'Version of the classifier to use',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snv', 'indel', 'sv', 'cnv'],
            doc => 'The type of variants being classified',
        }
    ],
    has => [
        prior_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'prior_result_id',
        },
    ],
    doc => '',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    unless($self->_validate_inputs) {
        my $err = $self->error_message;
        die $self->error_message('Failed to validate inputs: ' . $err);
    }

    unless($self->_prepare_staging_directory) {
        die $self->error_message('Failed to prepare staging directory.');
    }

    unless($self->_classify_variants) {
        die $self->error_message('Failed to run LOH.');
    }

    unless($self->_prepare_output_directory) {
        die $self->error_message('Failed to prepare output directory.');
    }

    unless($self->_promote_data) {
        die $self->error_message('Failed to promote data.');
    }

    $self->_reallocate_disk_allocation;

    return $self;
}

sub available_versions {
    return (1); #for yet unversioned things default to a version of 1
}

sub _validate_inputs {
    my $self = shift;

    unless($self->prior_result) {
        $self->error_message('No prior result found.');
        return;
    }

    my $path;
    if($self->prior_result->can('path')) {
        $path = $self->prior_result->path($self->variant_type . 's.hq.bed');
    } else {
        $path = join('/', $self->prior_result->output_dir, $self->variant_type . 's.hq.bed');
    }

    unless(-e $path) {
        $self->error_message('Could not find ' . $self->variant_type . ' file for prior result.');
        return;
    }

    my $version = $self->classifier_version;
    unless(grep($_ eq $version, $self->available_versions)) {
        $self->error_message('Unsupported classifier version passed.  Supported versions: ' . join(', ', $self->available_versions));
        return;
    }

    return 1;
}

1;
