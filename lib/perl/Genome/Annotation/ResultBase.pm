package Genome::Annotation::ResultBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::ResultBase {
    is_abstract => 1,
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        input_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
    has_param => [
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
    ],
};

sub output_filename {
    die "Abstract";
}

sub output_file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->output_filename);
}

sub input_result_file_path {
    my $self = shift;

    if ($self->input_result->can('output_file_path')) {
        return $self->input_result->output_file_path;
    } else {
        return $self->input_result->get_vcf($self->variant_type);
    }
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;
    $self->_run;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return File::Spec->join('/', 'model_data', 'software-result', $self->id);
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub param_names {
    my $self = shift;
    return ($self->is_many_param_names, $self->is_not_many_param_names);
}

sub is_many_param_names {
    my $self = shift;

    return map {$_->property_name} $self->__meta__->properties(is_many => 1, is_param => 1);
}

sub is_not_many_param_names {
    my $self = shift;

    my @names = map {$_->property_name} $self->__meta__->properties(is_many => 0, is_param => 1);
    my $names = Set::Scalar->new(@names);
    for my $name ($self->is_many_param_names, $self->is_many_input_names, $self->is_many_metric_names) {
        $names->delete($name . "_count");
        $names->delete($name . "_md5");
    }
    return $names->members();
}

sub param_hash {
    my $self = shift;

    my %hash;
    for my $param_name ($self->is_many_param_names) {
        $hash{$param_name} = [$self->$param_name];
    }
    for my $param_name ($self->is_not_many_param_names) {
        $hash{$param_name} = $self->$param_name;
    }
    return %hash;
}

sub metric_names {
    my $self = shift;

    return ($self->is_not_many_metric_names, $self->is_many_metric_names);
}

sub is_not_many_metric_names {
    my $self = shift;

    return map {$_->property_name} $self->__meta__->properties(is_many => 0, is_metric => 1);
}

sub is_many_metric_names {
    my $self = shift;

    return map {$_->property_name} $self->__meta__->properties(is_many => 1, is_metric => 1);
}

sub input_names {
    my $self = shift;
    return ($self->is_many_input_names, $self->is_not_many_input_names);
}

sub is_many_input_names {
    my $self = shift;

    return map {$_->property_name} $self->__meta__->properties(is_many => 1, is_input => 1);
}

sub is_not_many_input_names {
    my $self = shift;

    return map {$_->property_name} $self->__meta__->properties(is_many => 0, is_input => 1);
}

sub input_hash {
    my $self = shift;

    my %hash;
    for my $input_name ($self->is_many_input_names) {
        $hash{$input_name} = [$self->$input_name];
    }
    for my $input_name ($self->is_not_many_input_names) {
        $hash{$input_name} = $self->$input_name;
    }
    return %hash;
}
