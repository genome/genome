package Genome::Annotation::Detail::Result;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Detail::Result {
    is => 'Genome::SoftwareResult::Stageable',
};

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

    my @properties = $self->__meta__->properties(class_name => __PACKAGE__);
    return map {$_->property_name} grep {$_->is_param} @properties;
}

sub param_hash {
    my $self = shift;

    my %hash;
    for my $param_name ($self->param_names) {
        $hash{$param_name} = $self->$param_name;
    }
    return %hash;
}

