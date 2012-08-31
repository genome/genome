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

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class,
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

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
